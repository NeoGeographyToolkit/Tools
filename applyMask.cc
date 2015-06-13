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
#include <cmath>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/Cartography/GeoReference.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/Filter.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>

#include <asp/Core/Common.h>
#include <asp/Core/StereoSettings.h>

using namespace vw;

/// \file applyMask.cc Sets regions of an image as invalid according to a mask image.



/// Functor to evaluate if the sum squared value of a pixel is below a threshold
/// - In this program it is evaluating intersection error and generating a mask.
class InvalidPaintingFunctor  : public ReturnFixedType<PixelMask<PixelGray<float> > >
{
public:
  /// Construct with the thresholds.
  InvalidPaintingFunctor(const std::string &maskImagePath, 
                         const int imageWidth, const int imageHeight) : m_maskImage(maskImagePath)
  {
    // The mask image is loaded and we compute the scaling factors.
    int    maskHeight = m_maskImage.rows();
    int    maskWidth  = m_maskImage.cols();
    m_maskScaleX = static_cast<double>(maskWidth-1 ) / static_cast<double>(imageWidth);
    m_maskScaleY = static_cast<double>(maskHeight-1) / static_cast<double>(imageHeight);
    std::cout << "imageHeight = " << imageHeight << std::endl;
    std::cout << "imageWidth = " << imageWidth << std::endl;    
    std::cout << "maskHeight = " << maskHeight << std::endl;
    std::cout << "maskWidth = " << maskWidth << std::endl;
    std::cout << "scale x = " << m_maskScaleX << std::endl;
    std::cout << "scale y = " << m_maskScaleY << std::endl;
  }

  /// If the corresponding mask pixel is <= 0, make this pixel invalid.
  inline PixelMask<PixelGray<float> > operator() (PixelMask<PixelGray<float> > const &p, int32 col, int32 row, int32 plane=0) const
  {
  
    // Invalid pixels stay that way
    if (!is_valid(p))
      return p;  
  
    int maskImageCol = round(static_cast<double>(col) * m_maskScaleX);
    int maskImageRow = round(static_cast<double>(row) * m_maskScaleY);
    
    // Check for failure condition
    if ((maskImageCol >= m_maskImage.cols()) || (maskImageRow >= m_maskImage.rows()))
    {
      std::cout << "c = " << col << std::endl;
      std::cout << "r = " << row << std::endl;
      std::cout << "mC = " << maskImageCol << std::endl;
      std::cout << "mR = " << maskImageRow << std::endl;
      exit(0);
    }
    if (m_maskImage(maskImageCol, maskImageRow) > 0) // Indicates a valid pixel, keep it as is.
      return p;
    else // Change the pixel to invalid
    {
      PixelMask<PixelGray<float> > po(p);
      po.invalidate();
      return po;
    }
  }

private:

  DiskImageView<PixelGray<unsigned char> > m_maskImage;
  double m_maskScaleX;
  double m_maskScaleY;
};





template <class ImageT>
void write_geotiff_image( std::string                   const& filename,
                          vw::ImageViewBase<ImageT>     const& image,
                          vw::cartography::GeoReference const& georef,
                          double                        const  nodata_val,
                          asp::Options                  const& opt,
                          vw::TerminalProgressCallback  const& tpc ) {

  std::map<std::string, std::string> keywords;
  asp::block_write_gdal_image(filename, image.impl(), georef,
                              nodata_val, opt, tpc, keywords);
}

//TODO: Print help when no input arguments are used!
int main( int argc, char *argv[] ) {


  std::string inputImagePath, outputImagePath, maskImagePath;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h",        "Display this help message");

  asp::Options opt;
  general_options.add( asp::BaseOptionsDescription(opt) );

  po::options_description positional("");
  positional.add_options()
    ("input-image",   po::value(&inputImagePath),  "Path to input image file")
    ("mask-image",    po::value(&maskImagePath),   "Path to input mask file")
    ("output-image",  po::value(&outputImagePath), "Path to output image file");

  po::positional_options_description positional_desc;
  positional_desc.add("input-image",  1);
  positional_desc.add("mask-image",   1);
  positional_desc.add("output-image", 1);

  std::string usage("[options] <input-image> <mask-image> <output-image>\n");
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

  if ( !vm.count("input-image") || !vm.count("mask-image") || !vm.count("output-image"))
    vw_throw( vw::ArgumentErr() << "Requires input, mask, and output paths in order to proceed.\n\n"
              << usage << general_options );


  try {
  
    //TODO: Operate on multi-channel images of different data types!
    // Load the images from disk
    DiskImageView<PixelGray<float> > inputImage(inputImagePath);
  
  
    // Attempt to extract nodata value
    float noDataValue = -1;
    SrcImageResource *imageResource = DiskImageResource::open(inputImagePath);
    if ( imageResource->has_nodata_read() ) {
      noDataValue = imageResource->nodata_read();
      std::cout << "\t--> Extracted nodata value from file: " << noDataValue << ".\n";
    }
    else
      vw_throw( vw::ArgumentErr() << "Failed to extract NODATA value!\n\n");  
  

    // Copy georeference information from the input image
    vw::cartography::GeoReference georef;
    vw::cartography::read_georeference(georef, inputImagePath);
  
    // Set up thresholding function
    InvalidPaintingFunctor maskFunctor(maskImagePath, inputImage.cols(), inputImage.rows());
    write_geotiff_image(outputImagePath,
                        apply_mask( // Convert the mask back into pixel values
                                    vw::per_pixel_index_filter( // Add to mask using maskFunctor
                                                               create_mask( inputImage, // Set up mask
                                                                            noDataValue), 
                                                               maskFunctor),
                                    noDataValue),
                       georef, noDataValue, opt,
                       vw::TerminalProgressCallback( "applyMask", "Generating mask:")
                       );


/*
    // Set up output file
    boost::scoped_ptr<DiskImageResource> r(DiskImageResource::create(outputImagePath, outputView.format()));

    // Copy georeference information from the input image
    vw::cartography::GeoReference georef;
    vw::cartography::read_georeference(georef, inputImagePath);
    write_georeference( *r, georef );

    // write_georeferenced_image -> look up
    // Do everything!
    vw::block_write_image( *r, outputView,
                          vw::TerminalProgressCallback( "applyMask", "Generating mask:") );
*/
   
  } // End try statement
  
  
  catch (const Exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}










