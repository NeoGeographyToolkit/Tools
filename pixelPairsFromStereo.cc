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
#include <vw/InterestPoint.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>

// Not sure which of these we need but we get run-time errors otherwise...
#include <vw/InterestPoint.h>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <vw/Stereo/PreFilter.h>
#include <vw/Stereo/CorrelationView.h>
#include <vw/Stereo/CostFunctions.h>
#include <vw/Stereo/DisparityMap.h>
#include <asp/Core/DemDisparity.h>
#include <asp/Core/LocalHomography.h>

//#include <asp/Tools/stereo.h>
#include <stereo.h> // Using local copy

/// \file pixelPairsFromStereo.cc Generates a list of pixel point pairs from the output of the stereo tool





int main( int argc, char *argv[] ) {


  std::string inputImagePath, outputPath="";
  int pointSpacing=0;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h",        "Display this help message")  
    ("input-file,i",   po::value<std::string>(&inputImagePath),                "The input stereo file (xxx-D.tif)")
    ("output-file,o",  po::value<std::string>(&outputPath)->default_value(""), "Specify an output text file to store the program output")
    ("pointSpacing,p", po::value<int        >(&pointSpacing)->default_value(100), "Selected pixels are this far apart");
    
  po::positional_options_description positional_desc;

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options]" << std::endl << std::endl;
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
  
    std::ofstream outputFile(outputPath.c_str());
    if (outputFile.fail())
    {
      printf("Failed to open output file for writing\n");
      return 0;
    }
    
    printf("Loading image\n");

    vw::ImageViewRef<vw::PixelMask<vw::Vector2i> > inputImage = vw::DiskImageView<vw::PixelMask<vw::Vector2i> >(inputImagePath);

    printf("Done loading image\n");
    
    //if (!inputImage.is_valid_image())
    //{
    //  printf("Failed to load image!\n");
    //  return 0;
    //}
    printf("Image size: %d x %d\n", inputImage.rows(), inputImage.cols());

    vw::ImageViewRef<vw::PixelMask<vw::Vector2i> >::iterator iter = inputImage.begin();
    
    // First pass computes min, max, mean, and std_dev
    int count = 0;
    for (int row=0; row<inputImage.rows(); row+=1)
    {
    
      for (int col=0; col<inputImage.cols(); col+=1)
      {
        if ((row % pointSpacing == 0) && (col % pointSpacing == 0))
        {
          //printf("%d, %d\n", row, col);
          if (vw::is_valid(*iter)) // Skip invalid pixels
          {
            //printf("%d, %d\n", row, col);
            //std::cout << (*iter)[0] << std::endl;

            // Pixel in other image = original pixel location + stored offset
            int rightCol = floor(0.5 + col + (*iter)[0]);
            int rightRow = floor(0.5 + row + (*iter)[1]);

            // Write to file and record total
            outputFile << col <<","<< row <<","<< rightCol <<","<< rightRow << std::endl;
            ++count;
          } // End valid check
        } // End spacing check
        ++iter;
      } // End loop through cols

    } // End loop through rows

    // Done writing the output file
    outputFile.close();

    printf("Wrote %d pixel pairs\n", count);
    
  }
  catch (const vw::Exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}










