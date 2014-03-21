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
#include <asp/Core/InterestPointMatching.h>

#include <asp/IsisIO/IsisCameraModel.h>

//#include <asp/Tools/stereo.h>
#include <stereo.h> // Using local copy

/// \file stereoIpFind.h Generates a list of pixel pairs from two input images.





struct Parameters
{
  std::string leftFilePath;
  std::string rightFilePath;
  std::string outputPath;

  double elevationGuess;
};

bool parseInputParams(int argc, char *argv[], Parameters &params)
{

  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h",        "Display this help message")
    ("elevationGuess", po::value<double>(&params.elevationGuess)->default_value(0.0),
                       "Estimate of the elevation in the image");

  po::options_description positional("");
    positional.add_options()
      ("left",   po::value<std::string>(&params.leftFilePath),  "The left input file")
      ("right",  po::value<std::string>(&params.rightFilePath), "The right input file")
      ("output", po::value<std::string>(&params.outputPath),    "The output file");

    po::positional_options_description positional_desc;
    positional_desc.add("left",   1);
    positional_desc.add("right",  1);
    positional_desc.add("output", 1);

    std::string usage("[options] <left input path> <right input path> <output path>\n");
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

    if ( !vm.count("left") || !vm.count("right") || !vm.count("output") )
      vw_throw( vw::ArgumentErr() << "Requires <left input path> <right input path> <output path> input in order to proceed.\n\n"
                << usage << general_options );



}

bool produceInterestPoints(const Parameters &params)
{

  // Verify images are present
  boost::filesystem::path leftBoostPath (params.leftFilePath );
  boost::filesystem::path rightBoostPath(params.rightFilePath);

  if (!boost::filesystem::exists(boost::filesystem::path(params.leftFilePath)))
  {
    printf("Error: input file %s is missing!\n", params.leftFilePath.c_str());
    return false;
  }
  if (!boost::filesystem::exists(boost::filesystem::path(params.rightFilePath)))
  {
    printf("Error: input file %s is missing!\n", params.rightFilePath.c_str());
    return false;
  }


  // Load both images
  printf("Loading images left=%s and right=%s...\n",
         params.leftFilePath.c_str(),
         params.rightFilePath.c_str());
  vw::DiskImageView<vw::PixelGray<float> > leftDiskImage (params.leftFilePath );
  vw::DiskImageView<vw::PixelGray<float> > rightDiskImage(params.rightFilePath);

  //TODO: Error checks here?
  printf("Loading camera models...\n");
  vw::camera::IsisCameraModel leftCameraModel(params.leftFilePath);
  vw::camera::IsisCameraModel rightCameraModel(params.rightFilePath);

  //TODO: Where to get these from?
  double noData1 = -32767;
  double noData2 = -32767;

  // Set up Datum object for the moon, then apply the elevation guess
  vw::cartography::Datum datum("D_MOON");
  datum.set_semi_major_axis(datum.semi_major_axis() + params.elevationGuess);
  datum.set_semi_minor_axis(datum.semi_minor_axis() + params.elevationGuess);

  // Call function to find the matching interest points
  printf("Calling ipfind\n");
  bool ipFindResult = asp::ip_matching_w_alignment(&leftCameraModel, &rightCameraModel,
                                                   leftDiskImage,    rightDiskImage,
                                                   datum, params.outputPath,
                                                   noData1, noData2);

  return ipFindResult;
}


int main( int argc, char *argv[] )
{
  printf("Starting stereoIpFind\n");

  try {

    Parameters params;

    if (!parseInputParams(argc, argv, params))
      return 1;

    if (!produceInterestPoints(params))
      return 1;

    printf("stereoIpFind finished!\n");
  }
  catch (const vw::Exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}








