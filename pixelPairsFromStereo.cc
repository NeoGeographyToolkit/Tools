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

#include <vw/Math/Matrix.h>
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
  std::vector<float> leftAffineVals, rightAffineVals;

  // Affine transforms are identity matrix by default
  vw::Matrix<double> leftAffine, rightAffine;
  leftAffine  = vw::math::identity_matrix<3>();
  rightAffine = vw::math::identity_matrix<3>();
  // Note that the user inputs the matrices that were used in stereo,
  //  we need to compute the inverse of those matrices.

  po::options_description general_options("Options");
  general_options.add_options()
    ("help",        "Display this help message")
    ("pointSpacing", po::value<int        >(&pointSpacing)->default_value(100), "Selected pixels are this far apart")
    ("leftAffine",    po::value<std::vector<float> >(&leftAffineVals)->multitoken(),
                           "Affine transform which was applied to the left image pixels)")
    ("rightAffine",   po::value<std::vector<float> >(&rightAffineVals)->multitoken(),
                           "Affine transform which was applied to the right image pixels)");

  po::options_description positional("");
  positional.add_options()
    ("input-file",   po::value(&inputImagePath), "The input stereo file (xxx-D.tif)")
    ("output-file",  po::value(&outputPath),     "Specify an output text file to store the program output");

  po::positional_options_description positional_desc;
  positional_desc.add("input-file",  1);
  positional_desc.add("output-file", 1);

  std::string usage("[options] <input path> <output path>\n");
  po::variables_map vm;
  try {
    po::options_description all_options;
    all_options.add(general_options).add(positional);

    po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional_desc).style( po::command_line_style::unix_style ^ po::command_line_style::allow_short ).run(), vm );

    po::notify( vm );
  } catch (po::error const& e) {
    vw::vw_throw( vw::ArgumentErr() << "Error parsing input:\n"
                  << e.what() << "\n" << usage << general_options );
  }

  if ( !vm.count("input-file") || !vm.count("output-file") )
    vw_throw( vw::ArgumentErr() << "Requires <input path> and <output path> input in order to proceed.\n\n"
              << usage << general_options );

  const size_t NUM_AFFINE_ELEMENTS = 6; // The last three elements are 0 0 1
  if (vm.count("leftAffine"))
  {
    // Check input size
    if (leftAffineVals.size() != NUM_AFFINE_ELEMENTS)
      vw::vw_throw( vw::ArgumentErr() << "Incorrect number of left affine arguments passed in!\n"
                                      << usage << general_options );
    // Store affine values in a temp matrix object which is then inverted
    vw::Matrix<double> affineTemp = vw::math::identity_matrix<3>();
    for (int i=0; i<NUM_AFFINE_ELEMENTS; ++i)
      affineTemp.data()[i] = leftAffineVals[i];
    vw::vw_out() << "Read in left affine matrix:\n" << affineTemp << "\n";
    leftAffine = inverse(affineTemp);
  }
  if (vm.count("rightAffine"))
  {
    // Check input size
    if (rightAffineVals.size() != NUM_AFFINE_ELEMENTS)
      vw::vw_throw( vw::ArgumentErr() << "Incorrect number of right affine arguments passed in!\n"
                                      << usage << general_options );
    // Store affine values in a temp matrix object which is then inverted
    vw::Matrix<double> affineTemp = vw::math::identity_matrix<3>();
    for (int i=0; i<NUM_AFFINE_ELEMENTS; ++i)
      affineTemp.data()[i] = rightAffineVals[i];
    vw::vw_out() << "Read in right affine matrix:\n" << affineTemp << "\n";
    rightAffine = inverse(affineTemp);
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

            // Apply the affine transforms to the left and right pixels
            vw::Vector3 leftCoord (col,      row, 1);
            vw::Vector3 rightCoord(rightCol, rightRow, 1);
            vw::Vector3 leftTransformed  = leftAffine * leftCoord;
            vw::Vector3 rightTransformed = rightAffine * rightCoord;

            // Write to file and record total
            outputFile << leftTransformed [0] <<","<< leftTransformed [1] <<","
                       << rightTransformed[0] <<","<< rightTransformed[1] << std::endl;
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










