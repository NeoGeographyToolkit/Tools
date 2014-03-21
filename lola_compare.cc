// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/Cartography/GeoReference.h>
#include <vw/Math/Functors.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/Statistics.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>

using namespace vw;

/// \file lola_compare.cc Compares a DEM to a LOLA RDR file



// TODO: Move this class to a library file!
/// Class to compute a running standard deviation
/// - Code adapted from http://www.johndcook.com/standard_deviation.html
class RunningStatistics
{
public:
  RunningStatistics() : m_n(0), m_oldM(0), m_newM(0), m_oldS(0), m_newS(0) 
  {
    m_min = std::numeric_limits<double>::max();
    m_max = std::numeric_limits<double>::min();
  }

  void Clear()
  {
    m_n = 0;
  }

  void Push(double x)
  {
    m_n++;

    if (x < m_min)
        m_min = x;
    if (x > m_max)
        m_max = x;

    // See Knuth TAOCP vol 2, 3rd edition, page 232
    if (m_n == 1)
    {
      m_oldM = m_newM = x;
      m_oldS = 0.0;
    }
    else
    {
      m_newM = m_oldM + (x - m_oldM)/m_n;
      m_newS = m_oldS + (x - m_oldM)*(x - m_newM);

      // set up for next iteration
      m_oldM = m_newM; 
      m_oldS = m_newS;
    }
  }

  int NumDataValues() const
  {
    return m_n;
  }


  double Min() const
  {
    return m_min;
  }

  double Max() const
  {
    return m_max;
  }

  double Mean() const
  {
    return (m_n > 0) ? m_newM : 0.0;
  }

  double Variance() const
  {
    return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
  }

  double StandardDeviation() const
  {
    return sqrt( Variance() );
  }

private:
  int    m_n;
  double m_oldM, m_newM, m_oldS, m_newS;
  double m_max, m_min;

};

/// Function to write the output statistics to an output stream
bool writeOutput(const std::vector<float > &diffVector, 
                 const std::vector<double> &levels, 
                 const std::vector<size_t> &hist,
                 const RunningStatistics   &statCalc,
                       std::ostream        &stream)
{
  
  stream << "Mean elevation difference    = " << statCalc.Mean()              << std::endl;
  stream << "Standard deviation           = " << statCalc.StandardDeviation() << std::endl;
  stream << "Num intersecting LOLA points = " << statCalc.NumDataValues()     << std::endl << std::endl;

  // Print out the percentile distribution  
  stream << "Image distribution (5 percent intervals): " << std::endl;
  stream.setf(std::ios::fixed, std::ios::floatfield);
  stream.precision(2);
  stream.fill(' ');
  for (size_t i=0; i<hist.size(); ++i)
  {    
    float percent   = 0.05f * i;
    int   index      = diffVector.size() * percent;
    float percentile = diffVector[index];

    stream << "Percentile " << std::right << std::setw(4) << percent << " = " << std::setw(5) << percentile << std::endl;
  }
  stream << "Percentile " << std::right << std::setw(4) << std::left << 1.0f << " = " <<  diffVector.back() << std::endl;
  
  stream << std::endl; // Put a space between the "charts"

  // Print out the histogram
  stream << "Image histogram (20 bins):" << std::endl;
  stream.fill(' ');
  double ratio = 100.0 / statCalc.NumDataValues();
  for (size_t i=0; i<hist.size(); ++i)
  {    
    stream << std::right << std::setw(8) << levels[i] << " <-->" << std::setw(8) << levels[i+1] << " =" 
           << std::setw(6) << hist[i] << " =" << std::setw(5) << hist[i]*ratio << "%" << std::endl;
  }
  
  return true;
}

int main( int argc, char *argv[] )
{
  //This is pulled straight from the ASU DEM!  It may need to be different for other cases.
  const double MEAN_MOON_RADIUS = 1737400;

  const int numBins = 20; // 5% intervals

  std::string demPath, lolaDataPath, outputPath="", mapPath="", csvPath="";
  int removeHistogramOutliers=0;
  bool useAbsolute, sub180;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h",        "Display this help message")  
    ("output-file,o", po::value<std::string>(&outputPath)->default_value(""),     "Specify an output text file to store the program output")
    ("csv-log,c",     po::value<std::string>(&csvPath)->default_value(""),        "Specify an csv file to record individual point errors to")
    ("map,m",         po::value<std::string>(&mapPath)->default_value(""),        "Write output diagnostic image to file.")
    ("absolute",      po::bool_switch       (&useAbsolute)->default_value(false), "Output the absolute difference as opposed to just the difference.")
    ("sub180",        po::bool_switch       (&sub180)->default_value(false),      "Force all output degrees in csv file to be in -180 to 180 range.")
    ("limit-hist",    po::value<int        >(&removeHistogramOutliers)->default_value(0), "Limits the histogram to +/- N standard deviations from the mean");


  po::options_description positional("");
  positional.add_options()
    ("input-dem",   po::value(&demPath     ), "Path to DEM file")
    ("lola-points", po::value(&lolaDataPath), "Path to LOLA RDR point file");

  po::positional_options_description positional_desc;
  positional_desc.add("input-dem",   1);
  positional_desc.add("lola-points", 1);

  std::string usage("[options] <input DEM path> <input LOLA path>\n");
  po::variables_map vm;
  try {
    po::options_description all_options;
    all_options.add(general_options).add(positional);

    po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional_desc).style( po::command_line_style::unix_style ).run(), vm );

    po::notify( vm );
  } catch (po::error const& e) {
    vw::vw_throw( vw::ArgumentErr() << "Error parsing input:\n"
                  << e.what() << "\n" << usage << general_options );
  }

  if ( !vm.count("input-dem") || !vm.count("lola-points") )
    vw_throw( vw::ArgumentErr() << "Requires <input DEM path> and <input LOLA path> input in order to proceed.\n\n"
              << usage << general_options );



  // Load the input DEM
  ImageView<PixelGray<float> > inputDem;
  cartography::GeoReference georef;
  if (!read_georeferenced_image(inputDem, georef, demPath))
  {
    printf("Failed to read image!\n");
    return false;
  }
  std::cout << "Image size = " << inputDem.rows() << ", " << inputDem.cols() << std::endl;

  //// Pixel check
  //Vector2 pixelCoord(0,0);
  //Vector2 gdcCoord = georef.pixel_to_lonlat(pixelCoord);
  //std::cout << "0,0 = " << gdcCoord.x() << ", " << gdcCoord.y() << std::endl;

  //pixelCoord = Vector2(inputDem.cols()-1,inputDem.rows()-1);
  //gdcCoord = georef.pixel_to_lonlat(pixelCoord);
  //std::cout << inputDem.cols()-1 << "," << inputDem.rows()-1 << " = " << gdcCoord.x() << ", " << gdcCoord.y() << std::endl;

  // Loop through all lines in the LOLA data file
  std::ifstream lolaFile;
  lolaFile.open(lolaDataPath.c_str());
  std::string currentLine;
  std::getline(lolaFile, currentLine); // Skip header line
  if (lolaFile.fail())
  {
    printf("Failed to open LOLA data file %s\n", lolaDataPath.c_str());
    return false;
  }

  std::ofstream csvFile;
  if (csvPath.size() > 0) // Start CSV file and write header
  {
    csvFile.open(csvPath.c_str());
    csvFile << "# Latitude(Deg), Longitude(Deg), LOLA elevation(m), DEM elevation(m), difference(m)" << std::endl;
  }

  // Loop through all the lines in the file
  unsigned int numValidPixels = 0;
  RunningStatistics statCalc;
  double  demValue = 0;
  std::list<vw::Vector3> hitList;
  std::vector<float> diffVector;
  diffVector.reserve(5000); // Init this to an estimated size
  while (std::getline(lolaFile, currentLine))
  {
    // Get the location of the first few commas
    size_t timeComma    = currentLine.find(",");            
    size_t lonComma     = currentLine.find(",", timeComma+1); 
    size_t latComma     = currentLine.find(",", lonComma+1);  
    size_t radiusComma  = currentLine.find(",", latComma+1);  

    // Extract text containing the numbers we need    
    std::string lonString    = currentLine.substr(timeComma+1, lonComma   -timeComma-1);
    std::string latString    = currentLine.substr(lonComma +1, latComma   -lonComma -1);
    std::string radiusString = currentLine.substr(latComma +1, radiusComma-latComma -1);


    // Convert to floating point
    double ptLonDeg   = atof(lonString.c_str()   );
    double ptLatDeg   = atof(latString.c_str()   );
    double ptRadiusKm = atof(radiusString.c_str());

    // Convert to common elevation datum
    double lolaElevation = (ptRadiusKm*1000.0) - MEAN_MOON_RADIUS;

    // Find the equivalent location in the image
    Vector2 gdcCoord(ptLonDeg, ptLatDeg);
    Vector2 gccCoord   = georef.lonlat_to_point(gdcCoord);
    Vector2 pixelCoord = georef.point_to_pixel(gccCoord);

    //std::cout << "gdcCoord = " << gdcCoord << std::endl;
    //std::cout << "pixelCoord = " << pixelCoord << std::endl;

    // Throw out pixel misses
    if ( (pixelCoord[0] < 0) ||
         (pixelCoord[0] >= inputDem.cols()) ||
         (pixelCoord[1] < 0) ||
         (pixelCoord[1] >= inputDem.rows())   )
      continue; 

//    std::cout << "pixelCoord VALID = " << pixelCoord << std::endl;


    // Get DEM value and make sure it is valid
    demValue = inputDem(pixelCoord[0], pixelCoord[1])[0];
    if (demValue <= -32767) //TODO: Obtain the no-data value!
      continue; // Throw out flag values

//    std::cout << "demHeight     = " << demValue      << std::endl;
//    std::cout << "lolaElevation = " << lolaElevation << std::endl;    

    // Check difference
    float diff = static_cast<float>(demValue - lolaElevation);
    if (useAbsolute)
      diff = fabs(diff);
    ++numValidPixels;

    // Update records
    statCalc.Push(diff);
    diffVector.push_back(diff); // Record difference values for computing histogram later
    hitList.push_back(vw::Vector3(pixelCoord[0], pixelCoord[1], fabs(diff)));
    
    // Write a file listing all the compared locations
    if (csvPath.size() > 0)
    {
      csvFile.precision(12);
      if (sub180)
      {
      if (ptLonDeg > 180)
        ptLonDeg -= 360;
      }
      csvFile << ptLatDeg << ", " << ptLonDeg << ", " << lolaElevation << ", " << demValue << ", " << diff << std::endl;
    }
  } // End of loop through lines
  
  if (csvPath.size() > 0)
    csvFile.close();

  printf("Found %d valid pixels\n", numValidPixels);


  double meanValue = statCalc.Mean();
  double stdDevVal = statCalc.StandardDeviation();
  double minVal    = statCalc.Min();
  double maxVal    = statCalc.Max();

  if (removeHistogramOutliers > 0) // Cut off range at +/- N std
  {
    printf("Diff min = %lf\n", minVal);
    printf("Diff max = %lf\n", maxVal);

    minVal = meanValue - removeHistogramOutliers*stdDevVal;
    if (minVal < statCalc.Min())
      minVal = statCalc.Min();
    maxVal = meanValue + removeHistogramOutliers*stdDevVal;
    if (maxVal > statCalc.Max())
      maxVal = statCalc.Max();

    printf("Restricting histogram to +/- %d standard deviations: %lf <--> %lf\n", removeHistogramOutliers, minVal, maxVal);
  }

  double range   = maxVal - minVal;
  double binSize = range / numBins;
  double factor  = 1.0 / binSize;

  // Set up levels structure
  int numLevels = numBins + 1;
  std::vector<double> levels(numLevels);
  for (int i=0; i<numLevels; ++i)
    levels[i] = minVal + i*binSize;

  // Next fill in histogram
  std::vector<size_t> hist;
  hist.assign(numBins, 0);
  for (size_t i = 0; i < diffVector.size(); ++i)
  {
    float diff = diffVector[i];
    int   bin  = (int)floor( factor * (diff - minVal)  );
    if ((bin >= 0) && (bin < numBins)) // If removeHistogramOutliers is set some values will not fit in a bin
      ++(hist[bin]);

  } // End diff loop

  // Compute CDF distribution
  // - LOLA data is small enough to do this by sorting the values vector
  std::sort(diffVector.begin(), diffVector.end()); 

  // Write output to console
  writeOutput(diffVector, levels, hist, statCalc, std::cout);

  if (!outputPath.empty())  // Write output to file
  {
    std::ofstream file(outputPath.c_str());
    writeOutput(diffVector, levels, hist, statCalc, file);
    file.close();
  }

  // Optionally write a diagnostic image
  if (mapPath.size() > 0)
  {
    double outputFraction = 0.1; // Currently output image is just a fraction of the input image
    
    //TODO: Use existing functions to do this!
    // Initialize a diagnostic image
    // - Convert the input image from float32 grayscale to uint8 RGB
    int plotWidth  = inputDem.cols() * outputFraction;
    int plotHeight = inputDem.rows() * outputFraction;
    ImageView<PixelRGB<unsigned char> > diffPlot(plotWidth, plotHeight);
    
    for (int r=0; r<plotHeight; ++r)
    {
      for (int c=0; c<plotWidth; ++c)
      {
        //TODO: Get scaling factors for input image!
        float         originalValue = inputDem(r*outputFraction, c*outputFraction);
        unsigned char grayValue     = 0;//(originalValue - minVal) * (range/255.0);
        diffPlot(r,c) = PixelRGB<unsigned char>(grayValue, grayValue, grayValue);
      }
    }
    
    
    // Loop through all valid pixel hits and draw a scaled error in red
    std::list<vw::Vector3>::const_iterator iter;
    for (iter=hitList.begin(); iter!=hitList.end(); ++iter)
    {
      double diff = iter->z();
      unsigned char scaledDiff = (diff - minVal) * (range/255.0);
      diffPlot(iter->x(), iter->y()) = PixelRGB<unsigned char>(scaledDiff, 0, 0);
    }

    
    vw_out() << "Writing error diagram: " << mapPath << "\n";
    write_image(mapPath, diffPlot);

    
    
  } // End of wrriting diagnostic image

  return 0;
}






