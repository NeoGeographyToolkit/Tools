#include <iostream>
#include <fstream>
#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <asp/Core/Common.h>
#include <asp/Core/Macros.h>


/** 
  Simple tool to generate a point cloud from a DEM.
*/

using namespace std;
using namespace vw;

struct Options : public asp::BaseOptions {};

// Allows FileIO to correctly read/write these pixel types
namespace vw {
  typedef Vector<float64,6> Vector6;
  template<> struct PixelFormatID<Vector3>   { static const PixelFormatEnum value = VW_PIXEL_GENERIC_3_CHANNEL; };
  template<> struct PixelFormatID<Vector3f>  { static const PixelFormatEnum value = VW_PIXEL_GENERIC_3_CHANNEL; };
  template<> struct PixelFormatID<Vector4>   { static const PixelFormatEnum value = VW_PIXEL_GENERIC_4_CHANNEL; };
  template<> struct PixelFormatID<Vector6>   { static const PixelFormatEnum value = VW_PIXEL_GENERIC_6_CHANNEL; };
}

int main( int argc, char**argv ){

  if (argc < 3){
    std::cerr << "Usage: " << argv[0] << " input.tif output.tif" << std::endl;
    exit(1);
  }
  std::string input_name( argv[1] );
  std::string output_name( argv[2] );
  if ( input_name.empty() ) {
    std::cerr << "Please provide an input file." << std::endl;
    return 1;
  }

  cartography::GeoReference dem_georef;
  cartography::read_georeference( dem_georef, input_name );
  DiskImageView<float> dem( input_name );

  // To do: Make it bigger
  double nodata = -32768;
  boost::shared_ptr<DiskImageResource>
    dem_rsrc( new DiskImageResourceGDAL(input_name) );
  if (dem_rsrc->has_nodata_read()){
    nodata = dem_rsrc->nodata_read();
    cout<<"nodata =" << nodata << std::endl;
  }

  // Convert to cartesian
  ImageView<Vector3> point_cloud(dem.cols(), dem.rows());
  size_t count = 0;
  Vector3 mean_center;
  BBox2 bb;
  
  TerminalProgressCallback tpc("","");
  double inc_amount = 1.0 / double(dem.cols() );
  for ( int i = 0; i < dem.cols(); i++ ) {
    int local_count = 0;
    Vector3 local_mean;
      
    for (int j = 0; j < dem.rows(); j++ ) {

      if (dem(i, j) == nodata) continue;
      Vector2 lonlat = dem_georef.pixel_to_lonlat( Vector2(i,j) );
      bb.grow(lonlat);
      Vector3 lonlatrad( lonlat.x(), lonlat.y(), dem(i,j) );
      //std::cout.precision(20);
      //std::cout << "\nxyz " << i << ' ' << j << ' ' << lonlatrad << std::endl;
        
      Vector3 xyz = dem_georef.datum().geodetic_to_cartesian( lonlatrad );
      if ( xyz != Vector3() && xyz == xyz ) {
        point_cloud(i, j) = xyz;
        local_mean += xyz;
        local_count++;
      }
    }
    if ( local_count > 0 ) {
      local_mean /= double(local_count);
      double afraction = double(count) / double(count + local_count);
      double bfraction = double(local_count) / double(count + local_count);
      mean_center = afraction*mean_center + bfraction*local_mean;
      count += local_count;
    }
    tpc.report_incremental_progress( inc_amount );
  }
  tpc.report_finished();

  std::cout.precision(20);
  std::cout << "Found " << count << " valid points\n";
  std::cout << "Center is " << mean_center << std::endl;
  std::cout << "box is " << bb  << ' ' << bb.max() - bb.min() << std::endl;
  std::cout << "spacing is " << elem_quot(bb.max()-bb.min(), Vector2(point_cloud.cols()-1, point_cloud.rows()-1)) << std::endl;
  cout << "Writing: " << output_name << endl;
  Options opt;
  opt.gdal_options["BIGTIFF"] = "IF_SAFER";

  asp::block_write_gdal_image(output_name, point_cloud, opt,
                              TerminalProgressCallback("asp", "\t-->: "));


  return 0;
}
