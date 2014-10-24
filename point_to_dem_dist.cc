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

// TODO: Several of these functions are already present in pc_align

struct Options: asp::BaseOptions {
  // Input
  string csv, dem;
  // Output
  string output_prefix;
};

void handle_arguments( int argc, char *argv[], Options& opt ) {
  
  po::options_description general_options("");
  general_options.add_options()
    ("output-prefix,o", po::value(&opt.output_prefix), "Specify the output prefix.");

  general_options.add( asp::BaseOptionsDescription(opt) );

  po::options_description positional("");
  positional.add_options()
    ("csv", po::value(&opt.csv), "The csv file to find the distances from.")
    ("dem", po::value(&opt.dem), "The dem to the find distances to.");
  
  po::positional_options_description positional_desc;
  positional_desc.add("csv", 1);
  positional_desc.add("dem", 1);

  string usage("<csv> <dem>");
  bool allow_unregistered = false;
  std::vector<std::string> unregistered;
  po::variables_map vm =
    asp::check_command_line( argc, argv, opt, general_options, general_options,
                             positional, positional_desc, usage,
                             allow_unregistered, unregistered );

  if ( opt.dem.empty() || opt.csv.empty() )
    vw_throw( ArgumentErr() << "Missing input files.\n"
              << usage << general_options );

  // Set default output prefix if none was passed in
  if ( opt.output_prefix.empty() ) {
    opt.output_prefix = fs::basename(opt.csv) + "__" + fs::basename(opt.dem);
  }
  asp::create_out_dir(opt.output_prefix);

}

/// Helper function to check for parsing errors
void null_check(const char* token, string const& line){
  if (token == NULL)
    vw_throw( vw::IOErr() << "Failed to read line: " << line << "\n" );
}

/// Parse CSV file consisting of lat/lon/height triplets with an optional header.
void load_csv(string             const& file_name,
              cartography::Datum const& datum,
              vector<Vector3>         & lon_lat_height
              ){
  // Set up string copy buffer
  const int bufSize = 1024;
  char temp[bufSize];
  
  // Try to load the input file
  ifstream file( file_name.c_str() );
  if( !file ) {
    vw_throw( vw::IOErr() << "Unable to open file \"" << file_name << "\"" );
  }

  // Loop through all the lines of the file
  string line;
  char sep[] = ", \t";
  bool is_first_line = true;
  while ( getline(file, line, '\n') ){

    double lat, lon, height;
    
    strncpy(temp, line.c_str(), bufSize);
    const char* token = strtok(temp, sep); 
    null_check(token, line); // Throw if the token is null
    int ret = sscanf(token, "%lg", &lat);
    
    token = strtok(NULL, sep); 
    null_check(token, line);
    ret += sscanf(token, "%lg", &lon);
    
    token = strtok(NULL, sep); 
    null_check(token, line);
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

    lon_lat_height.push_back(Vector3(lon, lat, height));
  }

}

/// Compute the mean of an std::vector of doubles out to a length
double calc_mean(vector<double> const& errs, int len){
  if (len == 0) 
    return 0;
  double mean = 0.0;
  for (int i = 0; i < len; i++){
    mean += errs[i];
  }
  return mean/len;
}

typedef InterpolationView< EdgeExtensionView< ImageViewRef< PixelMask<float> >, 
                                        ConstantEdgeExtension>, 
                           BilinearInterpolation> InterpolationReadyDem;

/// Get ready to interpolate points on a DEM existing on disk.
InterpolationReadyDem load_interpolation_ready_dem(std::string           const& dem_path,
                                                   cartography::GeoReference  & georef) {
  // Load the georeference from the DEM
  bool is_good = cartography::read_georeference( georef, dem_path );
  if (!is_good) 
    vw_throw(ArgumentErr() << "DEM: " << dem_path << " does not have a georeference.\n");

  // Set up file handle to the DEM and read the nodata value
  DiskImageView<float> dem(dem_path);
  double nodata = numeric_limits<double>::quiet_NaN();
  boost::shared_ptr<DiskImageResource> dem_rsrc( new DiskImageResourceGDAL(dem_path) );
  if (dem_rsrc->has_nodata_read()) 
    nodata = dem_rsrc->nodata_read();

  // Set up interpolation + mask view of the DEM
  ImageViewRef< PixelMask<float> > masked_dem = create_mask(dem, nodata);
  return InterpolationReadyDem(interpolate(masked_dem));
}


/// Interpolates the DEM height at the input coordinate.
/// - Returns false if the coordinate falls outside the valid DEM area.
bool interp_dem_height(InterpolationReadyDem         const & dem, 
                       vw::cartography::GeoReference const & georef,
                       vw::Vector3                   const & lonlat,
                       double                              & dem_height) {
  // Convert the lon/lat location into a pixel in the DEM.
  vw::Vector2 pix = georef.lonlat_to_pixel(subvector(lonlat, 0, 2));
  double c = pix[0], 
         r = pix[1];
         
  // Quit if the pixel falls outside the DEM.
  if (c < 0 || c >= dem.cols()-1 || // TODO: This ought to be an image class function
      r < 0 || r >= dem.rows()-1 )
    return false;

  // Interpolate the DEM height at the pixel location
  vw::PixelMask<float> v = dem(c, r);
  if (!is_valid(v)) 
    return false;
    
  dem_height = v.child();
  return true;
}

/// As interp_dem_height but normalizes the longitude of input points.
bool interp_dem_height_safe(InterpolationReadyDem         const & dem, 
                            vw::cartography::GeoReference const & georef,
                            double                        const & mean_dem_lon,
                            vw::Vector3                   const & lonlat,
                            double                              & dem_height) {
  // Normalize the longitude to be within 360 degrees of the mean longitude
  vw::Vector3 lonlat_adjusted(lonlat);
  lonlat_adjusted[0] += 360.0*round((mean_dem_lon - lonlat[0])/360.0); // 360 deg adjustment
  
  return interp_dem_height(dem, georef, lonlat_adjusted, dem_height)
}

/// Computes the centroid lonlat point of a DEM
Vector2 compute_dem_mean_lonlat(vw::cartography::GeoReference const & georef,
                                int num_rows, int num_cols) {
  Vector2 mean_dem_lonlat =
    (georef.pixel_to_lonlat(Vector2(0, 0))
     + georef.pixel_to_lonlat(Vector2(num_cols-1, 0))
     + georef.pixel_to_lonlat(Vector2(0, num_rows-1)) 
     + georef.pixel_to_lonlat(Vector2(num_cols-1, num_rows-1))
     )/4.0;
  return mean_dem_lonlat;
}

//----------------------------------------------------------------------------

double calc_rmse(vector<double> const& errs, int len){
  double sum = 0.0;
  for (int i = 0; i < len; i++){
    sum += errs[i]*errs[i];
  }
  if (len == 0) return 0;
  return std::sqrt(sum/len);
}

int main( int argc, char *argv[] ){

  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    // Load the DEM and prepare it for interpolation
    cartography::GeoReference georef;
    InterpolationReadyDem masked_dem_interp(load_interpolation_ready_dem(opt.dem, georef));
    const int dem_num_rows = masked_dem_interp.rows();
    const int dem_num_cols = masked_dem_interp.cols();

    // Find the mean longitude for the DEM points. This will be used to adjust
    // the longitude of csv file points if need be.
    Vector2 mean_dem_lonlat = compute_dem_mean_lonlat(georef, dem_num_rows, dem_num_cols);
    double  mean_dem_lon    = mean_dem_lonlat[0];
      
    // Load the input CSV file into a vector of vectors
    vector<Vector3> lon_lat_height;
    load_csv(opt.csv, georef.datum(), lon_lat_height);
    std::cout << "Loaded: " << lon_lat_height.size() << " points from " << opt.csv << "." << std::endl;

    // Loop through all of the input points
    double nan = numeric_limits<double>::quiet_NaN();
    vector<double> valid_errors, all_dz(lon_lat_height.size(), nan), 
                   all_v(lon_lat_height.size(), nan), all_err(lon_lat_height.size(), nan);
    for (int i = 0; i < (int)lon_lat_height.size(); i++){

      Vector3 llh = lon_lat_height[i];
      double dem_height_here;
      if (!interp_dem_height_safe(masked_dem_interp, georef, mean_dem_lon, llh, dem_height_here))
        continue; // Skip if no DEM intersection
        
      // Compute and record the distance the input point height.
      double err = std::abs(llh[2] - dem_height_here);
      valid_errors.push_back(err);

      all_v[i] = v.child();
      all_dz[i] = dz;
      all_err[i] = err;
    }

    int len = valid_errors.size();
    vw_out() << "Computed " << len << " valid vertical point to dem errors." << std::endl;
    if (len == 0) 
      return 0;

    // Sort the error vector to make computing percentiles easy
    sort(valid_errors.begin(), valid_errors.end());

    // Compute and display statistics
    double p16 = valid_errors[std::min(len-1, (int)round(len*0.16))];
    double p50 = valid_errors[std::min(len-1, (int)round(len*0.50))];
    double p84 = valid_errors[std::min(len-1, (int)round(len*0.84))];
    vw_out() << "Error percentiles: "
             << " 16%: " << p16 << ", 50%: " << p50 << ", 84%: " << p84 << endl;
    double mean = calc_mean(valid_errors, len);
    vw_out() << "Mean error: " << mean << std::endl;
    double rmse = calc_rmse(valid_errors, len);
    vw_out() << "RMSE: " << rmse << std::endl;

    // Save the errors
    std::string output_file = opt.output_prefix + "-sample.csv";
    std::cout << "Writing: " << output_file << std::endl;
    ofstream outfile( output_file.c_str() );
    outfile << "# lat,lon,z_mHAE,z_samp_mHAE,dz_m,abs_dz_m" << endl;
    for (int i = 0; i < (int)lon_lat_height.size(); i++){
      Vector3 llh = lon_lat_height[i];
      outfile << std::setiosflags(ios::fixed) << std::setprecision(8) << llh[1] << ',' << llh[0] << ',' 
              << std::setiosflags(ios::fixed) << std::setprecision(3) << llh[2]
              << ',' << all_v[i] << ',' << all_dz[i] << ',' << all_err[i] << endl;
    }
    
  } ASP_STANDARD_CATCHES;
  
  return 0;
}

