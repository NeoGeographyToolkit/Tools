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

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <cstdlib>
#include <iostream>
#include <cmath>
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

struct Options : asp::BaseOptions {};

// Find the mean and std dev DEM of the given input sets of DEMs.
// Also output count.tif, showing for each pixel how many times
// it was encountered in the stack of DEMs.

int main( int argc, char *argv[] ){

  try {
    ImageView<double> mean_dem, std_dev_dem, count_dem;
    GeoReference georef;
    float nodata_val = numeric_limits<double>::quiet_NaN();

    if (argc <= 1){
      vw_throw(ArgumentErr() << "Missing input files.\n"
               << "Usage: " << argv[0] << " <dems>\n");
    }
  
    int start = 1;
    for (int i = start; i < argc; i++){
      std::string curr_file = argv[i];
      std::cout << "Reading: " << curr_file << std::endl;

      DiskImageResourceGDAL in_rsrc(curr_file);
      float curr_nodata_val = numeric_limits<double>::quiet_NaN();
      if ( in_rsrc.has_nodata_read() ) {
        curr_nodata_val = in_rsrc.nodata_read();
        vw_out() << "\tFound input nodata value: " << curr_nodata_val << std::endl;
      }else{
        vw_out() << "Nodata value not found in: " << curr_file << std::endl;
      }

      GeoReference curr_georef;
      read_georeference(curr_georef, in_rsrc);

      DiskImageView<float> curr_dem(in_rsrc);
      ImageViewRef< PixelMask<float> > masked_dem = create_mask(curr_dem, curr_nodata_val);

      InterpolationView<EdgeExtensionView< ImageViewRef < PixelMask<float> >, ConstantEdgeExtension>, BilinearInterpolation> masked_dem_interp = interpolate(masked_dem);


      if (i == start){
        mean_dem.set_size    ( curr_dem.cols(), curr_dem.rows() );
        std_dev_dem.set_size ( curr_dem.cols(), curr_dem.rows() );
        count_dem.set_size   ( curr_dem.cols(), curr_dem.rows() );
        georef = curr_georef;
        nodata_val = curr_nodata_val;
      }

      for (int col = 0; col < mean_dem.cols(); col++){
        for (int row = 0; row < mean_dem.rows(); row++){
          Vector2 pix(col, row);

          if (i != start){
            pix = curr_georef.lonlat_to_pixel(georef.pixel_to_lonlat(pix));
          }

          double x = pix[0], y = pix[1];
          double interp_val;
          bool is_in_box = ( 0 <= x && x <= masked_dem_interp.cols() - 1 &&
                             0 <= y && y <= masked_dem_interp.rows() - 1 );
          if (!is_in_box) continue;

          // The "if" statement below is to work around a bug in
          // interpolation with PixelMask.
          PixelMask<float> val;
          if (x == col && y == row)
            val = masked_dem(col, row);
          else
            val = masked_dem_interp(x, y);
          if (!is_valid(val)) continue;
          interp_val = val.child();

          mean_dem    (col, row) += interp_val;
          std_dev_dem (col, row) += interp_val*interp_val;
          count_dem   (col, row) += 1;
        }

      }
    }

    for (int col = 0; col < mean_dem.cols(); col++){
      for (int row = 0; row < mean_dem.rows(); row++){

        if (count_dem (col, row) == 0){
          mean_dem    (col, row) = nodata_val;
          std_dev_dem (col, row) = nodata_val;
          count_dem   (col, row) = nodata_val;
          continue;
        }

        double mean    = mean_dem(col, row)/count_dem(col, row);
        double std_dev = sqrt( std_dev_dem(col, row)/count_dem(col, row) - mean*mean );

        mean_dem    (col, row) = mean;
        std_dev_dem (col, row) = std_dev;
      }
    }

    Options opt;
    std::string mean_file = "mean.tif";
    std::cout << "Writing: " << mean_file << std::endl;
    asp::block_write_gdal_image(mean_file, pixel_cast<float>(mean_dem),
                                georef, nodata_val, opt);
  
    std::string std_dev_file = "std_dev.tif";
    std::cout << "Writing: " << std_dev_file << std::endl;
    asp::block_write_gdal_image(std_dev_file, pixel_cast<float>(std_dev_dem),
                                georef, nodata_val, opt);

    std::string count_file = "count.tif";
    std::cout << "Writing: " << count_file << std::endl;
    asp::block_write_gdal_image(count_file, pixel_cast<float>(count_dem),
                                georef, nodata_val, opt);

  } ASP_STANDARD_CATCHES;

  return 0;
}
