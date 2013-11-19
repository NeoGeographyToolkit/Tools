// __BEGIN_LICENSE__
//  Copyright (c) 2006-2012, United States Government as represented by the
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

// For every point in a given csv file draw a ray to the planet
// center, intersect the given dem (via interpolation), and compute
// the resulting distance from the starting point to the intersection
// point. Return the list of errors, and display statistics for them.

// The csv file is in the lat,lon,height format, with spaces or commas
// as separators. The first line may be a header file.

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/foreach.hpp>
namespace fs = boost::filesystem;
namespace po = boost::program_options;

#include <vw/Core/Functors.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/Statistics.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/tools/Common.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/Image/Interpolation.h>
#include <asp/Core/Macros.h>
#include <asp/Core/Common.h>

using namespace vw;
using namespace vw::cartography;
using namespace std;

struct Options: asp::BaseOptions {
  // Input
  string csv, dem;
  // Output
  //string output_prefix;
};

void handle_arguments( int argc, char *argv[], Options& opt ) {
  
  po::options_description general_options("");
  //general_options.add_options()
  //  ("output-prefix,o", po::value(&opt.output_prefix), "Specify the output prefix.");

  general_options.add( asp::BaseOptionsDescription(opt) );

  po::options_description positional("");
  positional.add_options()
    ("csv", po::value(&opt.csv),
     "The csv file to find the distances from.")
    ("dem", po::value(&opt.dem),
     "The dem to the find distances to.");
  
  po::positional_options_description positional_desc;
  positional_desc.add("csv", 1);
  positional_desc.add("dem", 1);

  string usage("<csv> <dem>");
  po::variables_map vm =
    asp::check_command_line( argc, argv, opt, general_options, general_options,
                             positional, positional_desc, usage );

  if ( opt.dem.empty() || opt.csv.empty() )
    vw_throw( ArgumentErr() << "Missing input files.\n"
              << usage << general_options );

  //if ( opt.output_prefix.empty() ) {
  //  opt.output_prefix =
  //    fs::basename(opt.csv) + "__" + fs::basename(opt.dem);
  // }
  //asp::create_out_dir(opt.output_prefix);

}

void null_check(const char* token, string const& line){
  if (token == NULL)
    vw_throw( vw::IOErr() << "Failed to read line: " << line << "\n" );
}

void load_csv(string const& file_name,
              cartography::Datum const& datum,
              vector<Vector3> & llh
              ){
  
  const int bufSize = 1024;
  char temp[bufSize];
  ifstream file( file_name.c_str() );
  if( !file ) {
    vw_throw( vw::IOErr() << "Unable to open file \"" << file_name << "\"" );
  }

  string line;
  char sep[] = ", \t";
  bool is_first_line = true;
  while ( getline(file, line, '\n') ){

    double lat, lon, height;
    
    strncpy(temp, line.c_str(), bufSize);
    const char* token = strtok(temp, sep); null_check(token, line);
    int ret = sscanf(token, "%lg", &lat);
    
    token = strtok(NULL, sep); null_check(token, line);
    ret += sscanf(token, "%lg", &lon);
    
    token = strtok(NULL, sep); null_check(token, line);
    ret += sscanf(token, "%lg", &height);
    
    // Be prepared for the fact that the first line may be the header.
    if (ret != 3){
      if (!is_first_line){
        vw_throw( vw::IOErr() << "Failed to read line: " << line << "\n" );
      }else{
        is_first_line = false;
        continue;
      }
    }
    is_first_line = false;

    llh.push_back(Vector3(lon, lat, height));
  }

}

double calc_mean(vector<double> const& errs, int len){
  double mean = 0.0;
  for (int i = 0; i < len; i++){
    mean += errs[i];
  }
  if (len == 0) return 0;
  return mean/len;
}

int main( int argc, char *argv[] ){

  Options opt;
  try {
    handle_arguments( argc, argv, opt );
    
    cartography::GeoReference georef;
    bool is_good = cartography::read_georeference( georef, opt.dem );
    if (!is_good) vw_throw(ArgumentErr() << "DEM: " << opt.dem
                           << " does not have a georeference.\n");

    DiskImageView<float> dem(opt.dem);
    double nodata = numeric_limits<double>::quiet_NaN();
    boost::shared_ptr<DiskImageResource> dem_rsrc
      ( new DiskImageResourceGDAL(opt.dem) );
    if (dem_rsrc->has_nodata_read()) nodata = dem_rsrc->nodata_read();

    ImageViewRef< PixelMask<float> > masked_dem = create_mask(dem, nodata);
    InterpolationView<EdgeExtensionView< ImageViewRef < PixelMask<float> >, ConstantEdgeExtension>, BilinearInterpolation> masked_dem_interp = interpolate(masked_dem);

    // Find the mean longitude for the DEM points. This will be used to adjust
    // the longitude of csv file points if need be.
    Vector2 mean_dem_lonlat =
      (georef.pixel_to_lonlat(Vector2(0, 0))
       + georef.pixel_to_lonlat(Vector2(dem.cols()-1, 0))
       + georef.pixel_to_lonlat(Vector2(0, dem.rows()-1)) 
       + georef.pixel_to_lonlat(Vector2(dem.cols()-1, dem.rows()-1))
       )/4.0;
    double mean_dem_lon = mean_dem_lonlat[0];
      
    vector<Vector3> llh;
    load_csv(opt.csv, georef.datum(), llh);
    std::cout << "Loaded: " << llh.size() << " points from " << opt.csv << "."
              << std::endl;

    vector<double> valid_errors;
    for (int i = 0; i < (int)llh.size(); i++){

      Vector3 p = llh[i];
      p[0] += 360.0*round((mean_dem_lon - p[0])/360.0); // 360 deg adjustment
      
      Vector2 pix = georef.lonlat_to_pixel(subvector(p, 0, 2));
      double c = pix[0], r = pix[1];
      if (c < 0 || c >= dem.cols()-1 || r < 0 ||
          r >= dem.rows()-1 || dem(c, r) == nodata) continue;

      PixelMask<float> v = masked_dem_interp(c, r);
      if (! is_valid(v)) continue;
      valid_errors.push_back(std::abs(p[2] - v.child()));
    }

    int len = valid_errors.size();
    vw_out() << "Computed " << len << " valid vertical point to dem errors."
             << std::endl;

    sort(valid_errors.begin(), valid_errors.end());
    if (len == 0) return 0;
    
    double p16 = valid_errors[std::min(len-1, (int)round(len*0.16))];
    double p50 = valid_errors[std::min(len-1, (int)round(len*0.50))];
    double p84 = valid_errors[std::min(len-1, (int)round(len*0.84))];
    vw_out() << "Error percentiles: "
             << " 16%: " << p16 << ", 50%: " << p50 << ", 84%: " << p84 << endl;
    
    double mean = calc_mean(valid_errors, len);
    vw_out() << "Mean error: " << mean << std::endl;
    
  } ASP_STANDARD_CATCHES;
  
  return 0;
}

