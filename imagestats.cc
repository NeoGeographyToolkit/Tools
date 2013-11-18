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
#include <vw/Image/ImageView.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/Statistics.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>


using namespace vw;

/// \file imagestats.cc Computes a number of statistics about an image



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
bool writeOutput(const std::vector<float > &cdfVector, 
                 const std::vector<double> &levels, 
                 const std::vector<size_t> &hist,
                 const RunningStatistics   &statCalc,
                       std::ostream        &stream)
{
  
  stream << "Mean elevation difference = " << statCalc.Mean()              << std::endl;
  stream << "Standard deviation        = " << statCalc.StandardDeviation() << std::endl;
  stream << "Num valid pixels          = " << statCalc.NumDataValues()     << std::endl << std::endl;

  // Print out the percentile distribution  
  stream << "Image distribution (approximated 5 percent intervals): " << std::endl;
  stream.setf(std::ios::fixed, std::ios::floatfield);
  stream.precision(2);
  stream.fill(' ');
  for (size_t i=0; i<cdfVector.size(); ++i)
  {    
    double percent = i * 0.05;
    stream << "Percentile " << std::right << std::setw(4) << percent << " = " << std::setw(7) << cdfVector[i] << std::endl;
  }
  
  stream << std::endl; // Put a space between the "charts"

  // Print out the histogram
  stream << "Image histogram (20 bins):" << std::endl;
  stream.fill(' ');
  double ratio = 100.0 / statCalc.NumDataValues();
  for (size_t i=0; i<hist.size(); ++i)
  {    
    stream << std::right << std::setw(8) << levels[i] << " <-->" << std::setw(7) << levels[i+1] << " =" 
           << std::setw(8) << hist[i] << " =" << std::setw(6) << hist[i]*ratio << "%" << std::endl;
  }
  
  return true;
}

//TODO: Print help when no input arguments are used!
int main( int argc, char *argv[] ) {

  const int numBins = 20; // 5% intervals

  std::string inputImagePath, outputPath="";
  int removeHistogramOutliers=0;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h",        "Display this help message")  
    ("input-file,i",  po::value<std::string>(&inputImagePath),                "Explicitly specify the input file")
    ("output-file,o", po::value<std::string>(&outputPath)->default_value(""), "Specify an output text file to store the program output")
    ("limit-hist",    po::value<int        >(&removeHistogramOutliers)->default_value(0), "Limits the histogram to +/- N standard deviations from the mean");
    
  po::positional_options_description positional_desc;
  positional_desc.add("input-image",  1);
  
  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <input-image>" << std::endl << std::endl;
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


  try {
  
    //TODO: Operate on multi-channel images of different data types!
    // Load the image from disk
    DiskImageView<PixelGray<float> > inputImage(inputImagePath);
  
    // First pass computes min, max, mean, and std_dev
    RunningStatistics statCalc;
    for (int row=0; row<inputImage.rows(); ++row)
    {
    
      for (int col=0; col<inputImage.cols(); ++col)
      {
      
        if (is_valid(inputImage(col,row))) // Skip invalid pixels
        {
          float diff = inputImage(col,row)[0];
          if (diff > -32767) // Avoid flag value
            statCalc.Push(diff);

//          if (diff > 500)
//          {
//            printf("Diff = %lf at row %d, col %d\n", diff, row, col);
//          }

        }
      } // End loop through cols

    } // End loop through rows
   
    double meanPixelValue   = statCalc.Mean();
    double stdDevPixelValue = statCalc.StandardDeviation();
    double minVal           = statCalc.Min();
    double maxVal           = statCalc.Max();

    if (removeHistogramOutliers > 0) // Cut off range at +/- N std
    {
      printf("Image min = %lf\n", minVal);
      printf("Image max = %lf\n", maxVal);

      minVal = meanPixelValue - removeHistogramOutliers*stdDevPixelValue;
      if (minVal < statCalc.Min())
        minVal = statCalc.Min();
      maxVal = meanPixelValue + removeHistogramOutliers*stdDevPixelValue;
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

    // CDF 
    vw::math::CDFAccumulator<float> cdfCalc(1000, 251); //TODO: What values to pass in?

    // Next pass fill in histogram
    std::vector<size_t> hist;
    hist.assign(numBins, 0);
    for (int row = 0; row < inputImage.rows(); row++)
    {
      for (int col = 0; col < inputImage.cols(); col++)
      {
        if (!is_valid(inputImage(col,row))) // Skip invalid pixels
          continue;
        float diff = inputImage(col,row)[0];
        if (diff <= -32767) // Avoid flag value
          continue;
        diff = fabs(diff);

        int bin = (int)floor( factor * (diff - minVal)  );
        if ((bin >= 0) && (bin < numBins)) // If removeHistogramOutliers is set some values will not fit in a bin
          ++(hist[bin]);
        cdfCalc(diff);

      } // End column loop
    } // End row loop

    // Fill out the CDF
    std::vector<float> cdfVector(numLevels);
    for (int i=0; i<numLevels; ++i) // --> 0.0 to 1.0
    {    
      double percent = 0.05 * i;
      cdfVector[i]   = cdfCalc.quantile(percent);
    }

    // Write output to console
    writeOutput(cdfVector, levels, hist, statCalc, std::cout);

    if (!outputPath.empty())  // Write output to file
    {
      std::ofstream file(outputPath.c_str());
      writeOutput(cdfVector, levels, hist, statCalc, file);
      file.close();
    }     
  }
  catch (const Exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}










