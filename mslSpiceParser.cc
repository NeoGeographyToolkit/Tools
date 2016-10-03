// __BEGIN_LICENSE__
//  Copyright (c) 2009-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file mslSpiceParser.cc
/** Tool for converting the MSL rover SPICE data into an easier to use format.

    TODO: Move this file into googlenasa repository?

    SPICE data is available here: http://naif.jpl.nasa.gov/pub/naif/pds/data/msl-m-spice-6-v1.0/mslsp_1000/data/
*/

#include <list>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <stdio.h>

// CSpice include files
#include "SpiceUsr.h"
#include "SpiceZfc.h"

#include <vw/Core/Exception.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Quaternion.h>


using namespace vw;
using std::endl;
using std::setprecision;
using std::setw;

struct Parameters
{
  std::vector<std::string>  kernelPaths;
  
 
  std::string traversePath; ///< Writes a nicely spaced CSV traverse trail 

  // Spacecraft Clock options
  std::string scQueryInputPath;  ///< List of times to query locations for
  std::string scQueryOutputPath; ///< output CSV path for time queries
  
  
  bool debug;
};

std::string usage = "Usage:  mslSpiceParser <SPICE folder>";


bool handle_arguments(int argc, char* argv[], Parameters &opt) 
{ 
  // TODO: Use boost program options to include some options!

  if (argc < 2)
  {
    printf("Error: Input directory is required!\n");
    std::cout << usage << std::endl;
    return false;
  }
  // Grab the input folder from the first argument
  //std::string prefix = "/home/smcmich1/data/mslSpice/";
  std::string prefix = argv[1];
  if (*prefix.rbegin() != '/')
    prefix += '/'; // Add trailing slash if needed
  
  if (argc > 2) // Set up for a list of Spacecraft Clock queries
  {
    opt.scQueryInputPath  = argv[2];
    opt.scQueryOutputPath = prefix + "queryLocations.csv";
  }
  else
    opt.scQueryInputPath = "";
  // TODO: Handle other args
  


  //opt.traversePath = "/home/smcmich1/data/mslSpice/spkTranslation.csv"; 
  opt.traversePath = prefix + "mslTraversePath.csv";

  // TODO: Is it safe to hard-code the following files???
  opt.debug = true;
  opt.kernelPaths.resize(8);
  opt.kernelPaths[0] = prefix + "spk/msl_surf_rover_tlm.bsp"; // This file has all the rover data to date
  opt.kernelPaths[1] = prefix + "fk/msl_v08.tf"; // Frame kernel
  opt.kernelPaths[2] = prefix + "lsk/naif0010.tls"; // Leap second kernel
  opt.kernelPaths[3] = prefix + "pck/pck00008.tpc"; // Planet kernel
  opt.kernelPaths[4] = prefix + "spk/msl_ls_ops120808_iau2000_v1.bsp"; // Rover sites relative to Mars
  opt.kernelPaths[5] = prefix + "spk/de425s.bsp";
  opt.kernelPaths[6] = prefix + "spk/mar085s.bsp";
  opt.kernelPaths[7] = prefix + "sclk/msl.tsc"; // Spacecraft clock


  return true;
}

bool dumpSpiceFile(const Parameters &params)  
{
  const int MARS_CODE = 499;
  const std::string J2000_FRAME_STRING  = "J2000";
  const std::string MARS_FRAME_STRING  = "IAU_MARS";
  //const std::string EARTH_STRING        = "earth";
  const std::string MARS_STRING        = "mars";

  const int MSL_CODE = -76;
  const int MSL_ROVER_CODE = -76000;
  // MSL related codes (such as instruments) are in th -76000 range
  // Note: MSL_SITE_1...399  ==>   -76501...-76899  (Fixed sites along its path)

  // Constants
  const double MEAN_MARS_RADIUS_KM = 3389.50; //TODO: Check this!
  const double RAD2DEG = 180/3.14159265359; //TODO!

  // Sample the position data at this interval
  const int DAYS_PER_INTERVAL = 90; // Convention of the MSL group
  const int STEPS_PER_DAY     = 24*10; // One position per ten minutes?
  const int NUM_STEPS         = DAYS_PER_INTERVAL * STEPS_PER_DAY; 

  // For the main traverse path, record a point every X meters to shorten the output CSV file
  const double MIN_RECORD_DIFF_METERS = 1;
  
  // TODO: Existing googlenasa stuff requires this output format:
  // sol#, lon, lat
  // -> Pick a good way to get the same "time" for each sol and record it.
  
  //TODO: Read these from file
  // Target and observer are common for all operations
  std::string target   = MARS_STRING;
  //std::string observer = EARTH_STRING;
  std::string absCorr  = "NONE";


/* TODO:

  - Get pointing data from the CK file?

*/

  // Load all of the kernel files
  printf("Loading all source kernels\n");
  size_t numKernels = params.kernelPaths.size();
  for (size_t i=0; i<numKernels; ++i)
  {
    furnsh_c(params.kernelPaths[i].c_str());
  }
  printf("Finished loading kernels.\n");


  SpiceChar timeString[51];
  SPICEDOUBLE_CELL(cover, 2000);
  SpiceDouble b, e;    
  

  // Set up output traverse file    
  std::ofstream outputFile;
  printf("Writing file %s\n", params.traversePath.c_str());
  outputFile.open(params.traversePath.c_str());
  outputFile << "# Longitude(deg), Latitude(deg), Elevation(m), Date" << std::endl;


  printf("\nSPK file coverage:\n");

  // Get list of bodies covered by the kernel
  SPICEINT_CELL(ids, 1000);
  std::string currentInputPath = params.kernelPaths[0]; // TODO: Clean this up!
  spkobj_c(currentInputPath.c_str(), &ids);

  // This should be a single value?

  // For each object -> In this case, the MSL and other fixed locations it has visited.
  for (int i=0; i<card_c(&ids); ++i)
  {
    // Get the current body ID and skip ahead until we hit the MSL list
    SpiceInt bodyId = SPICE_CELL_ELEM_I(&ids, i);

    // Get the coverage window
    scard_c(0, &cover);
    spkcov_c(currentInputPath.c_str(), bodyId, &cover);

    // Get the number of intervals in the coverage window.
    SpiceInt numIntervals = wncard_c(&cover);    

    //printf("Coverage for SPK object %ld -> %ld intervals\n", (long int)bodyId, (long int)numIntervals);

    if (bodyId != MSL_CODE) // We are only interested in the rover location
      continue;
    
    //return false;

    // Convert the coverage interval start and stop times to TDB calendar strings.
    SpiceDouble xyzVecLastRecorded[3];
    xyzVecLastRecorded[0] = 0;
    xyzVecLastRecorded[1] = 0;
    xyzVecLastRecorded[2] = 0;
    for (int j=0; j<numIntervals; j++)
    {
      // Get the endpoints of the jth interval.
      wnfetd_c(&cover, j, &b, &e);

      // Convert the endpoints to TDB calendar format time strings and display them.
      timout_c(b, "YYYY MON DD HR:MN:SC.### (TDB) ::TDB", 51, timeString);

      printf("\nInterval:  %ld\nStart:     %s\n", j, timeString);

      timout_c(e, "YYYY MON DD HR:MN:SC.### (TDB) ::TDB", 51, timeString);
      printf("Stop:      %s\n", timeString);

      //printf("b = %lf, e = %lf\n", b, e); 
      // Modify the position at even intervals
      SpiceDouble startEt = b+0; //// Hack to avoid different start times
      SpiceDouble stopEt  = e-0;
      SpiceDouble stepSize = (stopEt - startEt) / NUM_STEPS;

      // For the number of specified steps
      SpiceDouble lastEt;
      for (int i=0; i<NUM_STEPS; ++i)
      {

        SpiceDouble state[6];
        SpiceDouble lightTime;
        SpiceDouble et = startEt + stepSize*(SpiceDouble)i;
        lastEt = et;

        // Retrieve the position of the spacecraft at this time
        spkez_c(bodyId,  et,  MARS_FRAME_STRING.c_str(), absCorr.c_str(), MARS_CODE, state, &lightTime); 
        // Rover relative to Mars in J2000 frame, units are kilometers and km/sec
        if ( failed_c() )
        {
          printf("SpiceEditor Error: Failed to obtain SC J2000 position data at et %lf!!!!!!!\n", et);
          return false;
        }

        //SpiceDouble tol = 0; // Make sure we can do this 

        // Convert ephemeris time to spacecraft clock time
        //SpiceDouble sclkdp;
        //sce2c_c(LRO_CLOCK_ID, et, &sclkdp);

        // Convert ET to calendar time
        const SpiceInt STR_SIZE = 48;
        SpiceChar calendarString[STR_SIZE];
        etcal_c (et, STR_SIZE, calendarString);
        
        const int SPK_STATE_SIZE = 6; // X, Y, Z,  Vx, Vy, Vz
        //printf("State = ");
        for (int r=0; r<SPK_STATE_SIZE; ++r) {
           //printf("%lf    ", state[r]);
           //state[r] *= 1000.0; // Convert from km to m
        }
        //printf("\n");
         
        SpiceDouble radius, longitude, latitude;
        SpiceDouble xyzVec[3];
        vpack_c(state[0], state[1], state[2], xyzVec);
        reclat_c(xyzVec, &radius, &longitude, &latitude);
        double elevation = radius - MEAN_MARS_RADIUS_KM;
        //printf("---> %s # %lf, %lf, %lf\n", calendarString, longitude*RAD2DEG, latitude*RAD2DEG, elevation);
        

        // Compare the position to the previous position and reject if too close.
        double pointDist = sqrt( (state[0]-xyzVecLastRecorded[0])*(state[0]-xyzVecLastRecorded[0]) + 
                                 (state[1]-xyzVecLastRecorded[1])*(state[1]-xyzVecLastRecorded[1]) +
                                 (state[2]-xyzVecLastRecorded[2])*(state[2]-xyzVecLastRecorded[2])  );
        if (pointDist*1000 >= MIN_RECORD_DIFF_METERS) // Remember to convert km to meters
        {
          // Record the data to disk
          outputFile.precision(12);
          outputFile << longitude*RAD2DEG << ", " << latitude*RAD2DEG << ", " << elevation*1000 << ", " << calendarString << std::endl;
          xyzVecLastRecorded[0] = state[0];
          xyzVecLastRecorded[1] = state[1];
          xyzVecLastRecorded[2] = state[2];
        }
        
/*

        // Try to get the planet orientation at that time (planet_from_J2000)
        SpiceDouble  planet_from_J2000_R[3][3];
        pxform_c(J2000_FRAME_STRING.c_str(), MOON_FRAME_STRING.c_str(), et, planet_from_J2000_R);
        // This matrix converts J2000 orientations to LRO orientations at et
        if ( failed_c() )
        {
          printf("SpiceEditor Error: Failed to obtain planet J2000 orientation data at et %lf!!!!!!!!!!!!!!!!\n", et);
          return false;
        }

        SpiceDouble  newPosition[3];
        if (params.transformType == TRANSFORM_TYPE_LOCAL) // Apply existing lronac position offset
        {
          // Convert the LRONAC offset from meters to kilometers
          SpiceDouble lronacOffsetKm[3];
          lronacOffsetKm[0] = lronacOffset[0] / 1000.0;
          lronacOffsetKm[1] = lronacOffset[1] / 1000.0;
          lronacOffsetKm[2] = lronacOffset[2] / 1000.0;
          
          SpiceDouble instOffset_J2000Frame[3];
          // Now need to convert the LRONAC offset (spacecraft to frame) into J2000 coordinates
          // - Multiply by the inverted rotation matrix to go from spacecraft frame to J2000 frame
          mtxv_c(spacecraft_from_J2000_R, lronacOffsetKm, instOffset_J2000Frame);
          for (int r=0; r<3; ++r) // Add the rotated offset to the original coordinate
            newPosition[r] = state[r] + instOffset_J2000Frame[r]; // Adding km to km here
        }
        else // Apply global transform from file
        {
          // Convert the planet transform to the correct coordinate space
          SpiceDouble fixedPlanet_from_J2000_R[3][3];
          SpiceDouble fixedPlanet_in_J2000_T[3];

          // Convert correcting rotation and translation into J2000 frame from planet frame
          mxm_c(planetFixed_from_planet_R, planet_from_J2000_R, fixedPlanet_from_J2000_R);

          // Convert the state into units of meters
          SpiceDouble stateMeters[3]; // J2000 frame
          for (int r=0; r<3; ++r)
            stateMeters[r] = state[r]*1000;

          // Apply the corrective planet-based transform to the vector currently represented in J2000 space
          SpiceDouble stateMoon[3], stateMoonFixed[3], stateJ2000Fixed[3];
          mxv_c(planet_from_J2000_R, stateMeters, stateMoon);       // Convert to moon frame

          if (params.transformType == TRANSFORM_TYPE_GLOBAL)
          {
            // In the transform that comes out of lronacAngleSolver (courtesy of the AdjustedCamera class),
            //  the solved for rotation does not affect the position at all so just add the translation.
            // -  AdjustedCameraModel applies transform about camera, not about the planet!
            for (int r=0; r<3; ++r)
            {
              stateMoonFixed[r] = stateMoon[r] + planetFixed_from_planet_T[r];
            }
          }
          else // PC_ALIGN transform
          {
            // Apply the rotation to the position, then add in the translation
            mxv_c(planetFixed_from_planet_R, stateMoon, stateMoonFixed);
            for (int r=0; r<3; ++r)
            {
              stateMoonFixed[r] += planetFixed_from_planet_T[r];
            }
          }

          // Now that the position is corrected in body frame, convert back to J2000 frame
          mtxv_c(planet_from_J2000_R, stateMoonFixed, stateJ2000Fixed);
          for (int r=0; r<3; ++r) // Convert from meters back to kilometers
          {
            newPosition[r] = stateJ2000Fixed[r] / 1000.0;
          }

        } // End of global transform case
        */


      } // End of loop through steps

      outputFile.close(); // Close read data file


      //printf("Start time = %lf\n", startEt);
      //printf("Stop  time = %lf\n", lastEt);

    } // End loop through intervals for this body

    
  } // End of loop through bodies




  if (params.scQueryInputPath == "") // No queries passed in
    return true;


  // Done creating the traverse path, now handle specific time queries
  std::ifstream timeFile(params.scQueryInputPath.c_str());      
  printf("Writing file %s\n", params.scQueryOutputPath.c_str());
  outputFile.open(params.scQueryOutputPath.c_str());
  outputFile << "# Longitude(deg), Latitude(deg), Elevation(m), Name, scTime, Date" << std::endl;

  // For each object -> In this case, the MSL and other fixed locations it has visited.
  for (int i=0; i<card_c(&ids); ++i)
  {
    // Get the current body ID and skip ahead until we hit the MSL list
    SpiceInt bodyId = SPICE_CELL_ELEM_I(&ids, i);

    // Get the coverage window
    scard_c(0, &cover);
    spkcov_c(currentInputPath.c_str(), bodyId, &cover);

    // Get the number of intervals in the coverage window.
    SpiceInt numIntervals = wncard_c(&cover);    

    //printf("\nCoverage for SPK object %ld -> %ld intervals\n", (long int)bodyId, (long int)numIntervals);

    if (bodyId != MSL_CODE) // We are only interested in the rover location
      continue;
    
    //return false;

    // Convert the coverage interval start and stop times to TDB calendar strings.
    for (int j=0; j<numIntervals; j++)
    {
      // Get the endpoints of the jth interval.
      wnfetd_c(&cover, j, &b, &e);

      // Convert the endpoints to TDB calendar format time strings and display them.
      timout_c(b, "YYYY MON DD HR:MN:SC.### (TDB) ::TDB", 51, timeString);

      //printf("\nInterval:  %ld\nStart:     %s\n", j, timeString);

      timout_c(e, "YYYY MON DD HR:MN:SC.### (TDB) ::TDB", 51, timeString);
      //printf("Stop:      %s\n", timeString);

      //printf("b = %lf, e = %lf\n", b, e); 
      const SpiceDouble startEt = b;
      const SpiceDouble stopEt  = e;


      // Now loop through the times specified in the file
      std::string currentLine;
      while (std::getline(timeFile, currentLine))
      {
        //printf("%s\n", currentLine.c_str());
      
        // Treat #'s as comment lines
        if (currentLine[0] == '#')
          continue;
          
        // Extract the values from the line
        std::stringstream s(currentLine);
        std::string name;
        SpiceDouble spacecraftTime=0.0;
        s >> name >> spacecraftTime; //TODO: A more robust parsing of the text!
        name = name.substr(0, name.size()-1); // Strip trailing comma
        //printf("name = %s --- time = %lf\n", name.c_str(), spacecraftTime);

        // Convert from spacecraft ticks to ephemeris time
        SpiceDouble et = spacecraftTime; // It turns out the spacecraft clock time posted online is ET!
        //sct2e_c(bodyId, spacecraftTime, &et);
        //std::cout<< "SC time = " << spacecraftTime << " converted to ET = " << et << std::endl;

        if ((et < startEt) || (et > stopEt))
        {
          printf("WARNING: Input time %lf is outside the range spanned by the SPICE data!\n", spacecraftTime);
          outputFile << 0 << ", " << 0 << ", " << 0 << ", \"" << name << "\", " << spacecraftTime << ", \"\"" << std::endl; // Write a blank row
          continue;
        }

        // Retrieve the position of the spacecraft at this time
        SpiceDouble state[6];
        SpiceDouble lightTime;
        spkez_c(bodyId,  et,  MARS_FRAME_STRING.c_str(), absCorr.c_str(), MARS_CODE, state, &lightTime); 
        // Rover relative to Mars in J2000 frame, units are kilometers and km/sec
        if ( failed_c() )
        {
          printf("SpiceEditor Error: Failed to obtain SC J2000 position data at et %lf!!!!!!!\n", et);
          outputFile << 0 << ", " << 0 << ", " << 0 << ", \"" << name << "\", " << spacecraftTime << ", \"\"" << std::endl; // Write a blank row
          continue;
        }

        // Convert ET to calendar time
        const SpiceInt STR_SIZE = 48;
        SpiceChar calendarString[STR_SIZE];
        etcal_c (et, STR_SIZE, calendarString);
        
        const int SPK_STATE_SIZE = 6; // X, Y, Z,  Vx, Vy, Vz
        //printf("State = ");
        for (int r=0; r<SPK_STATE_SIZE; ++r) {
           //printf("%lf    ", state[r]);
           //state[r] *= 1000.0; // Convert from km to m
        }
        //printf("\n");
         
        SpiceDouble radius, longitude, latitude;
        SpiceDouble xyzVec[3];
        vpack_c(state[0], state[1], state[2], xyzVec);
        reclat_c(xyzVec, &radius, &longitude, &latitude);
        double elevation = radius - MEAN_MARS_RADIUS_KM;
        //printf("---> %s # %lf, %lf, %lf\n", calendarString, longitude*RAD2DEG, latitude*RAD2DEG, elevation);
        
        // Record the data to disk
        outputFile.precision(12);
        outputFile << longitude*RAD2DEG << ", " << latitude*RAD2DEG << ", " << elevation*1000 << ", \"" 
                   << name << "\", " << spacecraftTime << ", \"" << calendarString << "\"" << std::endl;
        
      } // End of loop through requested SC times

    } // End loop through intervals for this body
    
  } // End of loop through bodies

  outputFile.close(); // Close data files
  timeFile.close();


  return true;
}




int main(int argc, char* argv[]) 
{
  // Parse the input parameters
  Parameters params;
  if (!handle_arguments(argc, argv, params))
  {
    printf("Failed to parse input parameters!\n");
    return false;
  }

  dumpSpiceFile(params);

  printf("Exiting spkDump program.\n");
  return 0;
}






