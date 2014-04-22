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

struct Options : asp::BaseOptions {};

// Convolve a DEM with exp(-sigma*x^2). The input DEM must
// have its invalid pixels masked.

template <class ImageT>
class BlurDEM:
  public ImageViewBase< BlurDEM<ImageT> > {
  ImageT m_img;
  int m_search_dist; // half of size of window to convolve with
  ImageView<double> m_gauss_kernel;
  typedef typename ImageT::pixel_type PixelT;

public:

  typedef PixelT pixel_type;
  typedef PixelT result_type;
  typedef ProceduralPixelAccessor<BlurDEM> pixel_accessor;

  BlurDEM( ImageViewBase<ImageT> const& img,
                     double blur_sigma) :
    m_img(img.impl()) {
    VW_ASSERT(blur_sigma > 0,
              ArgumentErr() << "Expecting positive sigma.");

    // Cut the gaussian exp(-sigma*x^2) where its value is 'scale'.
    double scale = 0.001;
    m_search_dist = (int)ceil(sqrt(-log(scale)/blur_sigma));
    std::cout << "Search distance is " << m_search_dist << std::endl;

    // The gaussian kernel
    int h = m_search_dist;
    m_gauss_kernel.set_size(2*h+1, 2*h+1);
    for (int c = 0; c < m_gauss_kernel.cols(); c++){
      for (int r = 0; r < m_gauss_kernel.rows(); r++){
        double r2 = double(c-h)*(c-h) + double(r-h)*(r-h);
        m_gauss_kernel(c, r) = exp(-blur_sigma*r2);
      }
    }
  
  }

  inline int32 cols() const { return m_img.cols(); }
  inline int32 rows() const { return m_img.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor(*this); }

  inline result_type operator()( size_t i, size_t j, size_t p=0 ) const {
    vw_throw( NoImplErr() << "BlurDEM: operator() not implemented.\n" );
  }
  
  
  typedef CropView< ImageView<PixelT> > prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {

    // Crop into an expanded box as to have enough pixels to do
    // the blurring at every pixel in the current box.
    int h = m_search_dist; // shorten
    BBox2i biased_box = bbox;
    biased_box.expand(h+1);
    biased_box.crop(bounding_box(m_img));
    ImageView<PixelT> img( crop( m_img, biased_box ) );
    ImageView<PixelT> filled_img = copy(img);

    int nc = img.cols(), nr = img.rows(); // shorten
    for (int row = 0; row < nr; row++){
      for (int col = 0; col < nc; col++){
        PixelT V; V.validate();
        double sum = 0.0;
        for (int c = std::max(col-h, 0); c <= std::min(col+h, nc-1); c++){
          for (int r = std::max(row-h, 0); r <= std::min(row+h, nr-1); r++){
            if (!is_valid(img(c, r))) continue;
            double wt = m_gauss_kernel(c-col+h, r-row+h);
              V   += wt*img(c, r);
              sum += wt;
          }
        }
        if (sum > 0) filled_img(col, row) = V/sum;
        
      }
    } 
      
    return prerasterize_type(filled_img,
                             -biased_box.min().x(), -biased_box.min().y(),
                             cols(), rows());
      
  }
  template <class ImgT>
  inline void rasterize( ImgT const& img, BBox2i const& bbox ) const {
    vw::rasterize( prerasterize(bbox), img, bbox );
  }
  
};

template <class ImgT>
BlurDEM<ImgT>
blur_dem( ImageViewBase<ImgT> const& img, double blur_sigma) {
  typedef BlurDEM<ImgT> result_type;
  return result_type( img.impl(), blur_sigma);
}

int main( int argc, char *argv[] ){

  Options opt;

  if (argc < 4){
    std::cerr << "Usage: " << argv[0] << " input-DEM.tif blur_sigma output-DEM.tif"
              << std::endl;
    exit(1);
  }
  
  std::string infile = argv[1];
  double blur_sigma = atof(argv[2]);
  std::string outfile = argv[3];
  
  std::cout << "Reading: " << infile << std::endl;
  std::cout << "blur sigma is " << blur_sigma << std::endl;
  
  DiskImageResourceGDAL in_rsrc(infile);
  float nodata_val = -32768;
  if ( in_rsrc.has_nodata_read() ) {
    nodata_val = in_rsrc.nodata_read();
    vw_out() << "\tFound input nodata value: " << nodata_val << std::endl;
  }else{
    std::cerr << "Nodata value not found in: " << infile << std::endl;
    exit(1);
  }
  
  DiskImageView<float> dem(in_rsrc);


  GeoReference georef;
  read_georeference(georef, in_rsrc);
  
  std::cout << "Writing: " << outfile << std::endl;
  block_write_gdal_image(outfile,
                         apply_mask
                         (blur_dem
                          (create_mask(dem, nodata_val), blur_sigma),
                          nodata_val),
                         georef, nodata_val, opt,
                         TerminalProgressCallback("asp","")
                         );
  return 0;
  
}
