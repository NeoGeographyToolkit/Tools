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

  std::string demPath, lolaDataPath, outputPath="";
  int removeHistogramOutliers=0;
  bool useAbsolute;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h",        "Display this help message")  
    ("input-dem,d",   po::value<std::string>(&demPath),                           "Explicitly specify the DEM  file")
    ("input-lola,l",  po::value<std::string>(&lolaDataPath),                      "Explicitly specify the LOLA file")
    ("output-file,o", po::value<std::string>(&outputPath)->default_value(""),     "Specify an output text file to store the program output")
    ("absolute",      po::value<bool       >(&useAbsolute)->default_value(false), "Output the absolute difference as opposed to just the difference.")
    ("limit-hist",    po::value<int        >(&removeHistogramOutliers)->default_value(0), "Limits the histogram to +/- N standard deviations from the mean");
    
  po::positional_options_description positional_desc;
  positional_desc.add("input-image", 1);
  positional_desc.add("input-lola",  1);
  
  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <input-DEM> <input-LOLA>" << std::endl << std::endl;
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(general_options).positional(positional_desc).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << usage.str();
    return 1;
  }

  ImageView<PixelGray<float> > inputDem;
  cartography::GeoReference georef;

  if (!read_georeferenced_image(inputDem, georef, demPath))
  {
    printf("Failed to read image!\n");
    return false;
  }
  std::cout << "Image size = " << inputDem.rows() << ", " << inputDem.cols() << std::endl;
    
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
  
  // Loop through all the lines in the file
  unsigned int numValidPixels = 0;
  RunningStatistics statCalc;
  double  demValue = 0;
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
    
    if ( (pixelCoord[0] < 0) ||
         (pixelCoord[0] >= inputDem.cols()) ||
         (pixelCoord[1] < 0) ||
         (pixelCoord[1] >= inputDem.rows())   )
      continue; // Throw out pixel misses

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

  } // End of loop through lines

  printf("Found %d valid pixels\n", numValidPixels);


  double meanValue = statCalc.Mean();
  double stdDevVal = statCalc.StandardDeviation();
  double minVal    = statCalc.Min();
  double maxVal    = statCalc.Max();

  if (removeHistogramOutliers > 0) // Cut off range at +/- N std
  {
    printf("Image min = %lf\n", minVal);
    printf("Image max = %lf\n", maxVal);

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
  
  return 0;
}






