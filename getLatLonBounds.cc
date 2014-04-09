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

/// \file getLatLonBounds.cc Estimates the latitude and longitude bounding box.




int main( int argc, char *argv[] )
{
  std::string inputImagePath;
//  double elevationGuess;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h", "Display this help message")  
    /*("elevationGuess", po::value<double>(&elevationGuess)->default_value(0),  "Specify an elevation to use when computing the coordinates")*/;


  po::options_description positional("");
  positional.add_options()
    ("inputImage",   po::value(&inputImagePath), "Path to input geotiff file");

  po::positional_options_description positional_desc;
  positional_desc.add("inputImage", 1);
  
  std::string usage("[options] <input image>\n");
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

  if ( !vm.count("inputImage") )
    vw_throw( vw::ArgumentErr() << "Requires <input image> in order to proceed.\n\n"
              << usage << general_options );



  // Load the input DEM georeference
  cartography::GeoReference georef; 
  boost::scoped_ptr<DiskImageResource> inputImage(DiskImageResource::open(inputImagePath));
  if (!read_georeference(georef, *inputImage))
  {
    //vw_out() << "Failed to read input image!\n";
    std::cout << "Failed to read input image georeference!\n";
    return false;
  }

  // Compute the GDC bounding box of the entire image
  BBox2 fullImageBox(Vector2(0,0), Vector2(inputImage->cols(), inputImage->rows()));
  BBox2 gdcBox = georef.pixel_to_lonlat_bbox(fullImageBox);

  // Print out the results
  /*vw_out() <<   "Min latitude  = " << gdcBox.min()[1]
           << "\nMax latitude  = " << gdcBox.max()[1]
           << "\nMin longitude = " << gdcBox.min()[0]
           << "\nMax longitude = " << gdcBox.max()[0];
 */

  std::cout <<   "Min latitude  = " << gdcBox.min()[1]
            << "\nMax latitude  = " << gdcBox.max()[1]
            << "\nMin longitude = " << gdcBox.min()[0]
            << "\nMax longitude = " << gdcBox.max()[0] << "\n";


  return 0;
}






