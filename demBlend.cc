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

/// \file demBlend.cc
///

// A tool to mosaic and blend DEMs, and output the mosaic as tiles.
// Unless USE_GRASSFIRE is set to 0 below, it expects DEMs to have an
// alpha channel with the grassfire weights, which are used for
// blending the DEMs. Such a DEM can be obtained from a regular DEM
// using the VisionWorkbench grassfirealpha command.

// Note 1: In practice, the tool may be more efficient if the entire
// mosaic is written out as one single large image, rather than being
// broken up into tiles. To achieve that, just specify to the tool a
// very large tile size, and use 0 for the tile index in the command
// line options.

// Note 2: The tool can be high on memory usage, so processes for
// individual tiles may need to be run on separate machines.

// To do:
// Deal with the protobuf dependency in the build system.
// Fix cmake to build in build/ and not in base dir.
// Add unit tests.

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <limits>
using namespace std;

#include <vw/FileIO.h>
#include <vw/Image.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;

// TODO: Why is this not installed by ASP?
//#include <asp/Core/AntiAliasing.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;



// TODO: Move all of these functions to vision workbench!

// Copy-pasted from PhotomotryTK ------

//TODO: Replace with georef BBox function
Vector4 ComputeGeoBoundary(cartography::GeoReference Geo, int width, int height){

  // Get the lonlat coordinates of the four pixels corners of the image.

  Vector4 corners;
  Vector2 leftTopPixel(0,0);
  Vector2 leftTopLonLat = Geo.pixel_to_lonlat(leftTopPixel);

  Vector2 rightBottomPixel(width-1, height-1);
  Vector2 rightBottomLonLat = Geo.pixel_to_lonlat(rightBottomPixel);

  float minLon = leftTopLonLat(0);
  float minLat = leftTopLonLat(1);
  float maxLon = rightBottomLonLat(0);
  float maxLat = rightBottomLonLat(1);

  if (maxLat<minLat){
    float temp = minLat;
    minLat = maxLat;
    maxLat = temp;
  }

  if (maxLon<minLon){
    float temp = minLon;
    minLon = maxLon;
    maxLon = temp;
  }

  corners(0) = minLon;
  corners(1) = maxLon;
  corners(2) = minLat;
  corners(3) = maxLat;

  return corners;
}

//TODO: replace
vw::Vector4 getImageCorners(std::string imageFile){

  // Get the four corners of an image, that is the lon-lat coordinates of
  // the pixels in the image corners.

  // Note: Below we assume that the image is uint8. In fact, for the
  // purpose of calculation of corners the type of the image being
  // read does not matter.
  DiskImageView<PixelMask<PixelGray<uint8> > >  image(imageFile);

  cartography::GeoReference imageGeo;
  bool is_good = read_georeference(imageGeo, imageFile);
  if (!is_good){
    std::cerr << "No georeference found in " << imageFile << std::endl;
    exit(1);
  }
  Vector4 imageCorners = ComputeGeoBoundary(imageGeo, image.cols(), image.rows());
  return imageCorners;
}

// If prefix is "dir/out", create directory "dir"
void create_out_dir(std::string out_prefix){

  boost::filesystem::path out_prefix_path(out_prefix);
  if (out_prefix_path.has_parent_path()) {
    if (!boost::filesystem::is_directory(out_prefix_path.parent_path())) {
      vw_out() << "\nCreating output directory: "
               << out_prefix_path.parent_path() << std::endl;
      boost::filesystem::create_directory(out_prefix_path.parent_path());
    }
  }

  return;
}

// ------

// TODO: This is needed because the existing subsample view does not
//        divide odd sizes as would be expected.
// TODO: Why doesn't this work?
// Simple 2x downsample class
template <class ImageT>
class Downsample2xView : public ImageViewBase<Downsample2xView<ImageT> > {
  ImageT m_child;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename boost::remove_reference<typename ImageT::result_type>::type result_type;
  typedef ProceduralPixelAccessor<Downsample2xView<ImageT> > pixel_accessor;

  Downsample2xView( ImageT const& image ) : m_child(image) {}

  inline int32 cols  () const { return m_child.cols()/2; }
  inline int32 rows  () const { return m_child.rows()/2; }
  inline int32 planes() const { return m_child.planes(); }
  inline pixel_accessor origin() const { return pixel_accessor( *this ); }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const {
    return m_child(i*2,j*2,p);
  }

  typedef Downsample2xView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    return prerasterize_type(m_child.prerasterize(bbox)); }
  template <class DestT>
  inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
    vw::rasterize( prerasterize(bbox), dest, bbox );
  }
};

template <class ViewT>
Downsample2xView<ViewT>
downsample2x( ImageViewBase<ViewT> const& view ) {
  return Downsample2xView<ViewT>( view.impl() );
}

// Simple 2x upsample class
template <class ImageT>
class Upsample2xView : public ImageViewBase<Upsample2xView<ImageT> > {
  ImageT m_child;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename boost::remove_reference<typename ImageT::result_type>::type result_type;
  typedef ProceduralPixelAccessor<Upsample2xView<ImageT> > pixel_accessor;

  Upsample2xView( ImageT const& image ) : m_child(image) {}

  inline int32 cols  () const { return m_child.cols()*2; }
  inline int32 rows  () const { return m_child.rows()*2; }
  inline int32 planes() const { return m_child.planes(); }
  inline pixel_accessor origin() const { return pixel_accessor( *this ); }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const {
    if ((i % 2 != 0) || (j % 2 != 0)) // Return zero except on even indices
      return result_type();
    // Otherwise compute image location
    return m_child(i/2,j/2,p);
  }

  typedef Upsample2xView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    return prerasterize_type(m_child.prerasterize(bbox)); }
  template <class DestT>
  inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
    vw::rasterize( prerasterize(bbox), dest, bbox );
  }
};

template <class ViewT>
Upsample2xView<ViewT>
upsample2x( ImageViewBase<ViewT> const& view ) {
  return Upsample2xView<ViewT>( view.impl() );
}

//TODO: Speed optimization
/// Generate a laplacian image pyramid for the given image
template <typename T>
bool buildLaplacianImagePyramid(const ImageViewBase<T> &inputImage,
                                const int numLevels,
                                std::vector<ImageView<typename T::pixel_type> > &pyramid)
{
  // Set up the gaussian smoothing kernel
  // - Copied from stereo PyramidCorrelationView, experiment with this.
  // - This is a 5x5 seperable kernel.
  std::vector<typename DefaultKernelT<typename T::pixel_type>::type > kernel(5);
  kernel[0] = kernel[4] = 1.0/16.0;
  kernel[1] = kernel[3] = 4.0/16.0;
  kernel[2] = 6.0/16.0;

  // Initialize output pyramid levels
  pyramid.resize(numLevels);

  // Initialize temporary storage for blurred images
  // TODO: Any way to avoid allocating these buffers?
  std::vector<ImageView<typename T::pixel_type> > blurredImages(numLevels);
  
  // First create a full Gaussian pyramid
  // Loop through pyramid levels
  printf("Building Gaussian pyramid\n");
  blurredImages[0] = inputImage;
  
  std::stringstream path;
  path << "/home/smcmich1/data/demBlendTest/pyrBlur" << 0 << ".tiff";
  write_image(path.str(), blurredImages[0]);  
  
  //std::cout << "blurred: " << blurredImages[0].cols() << " " << blurredImages[0].rows() << std::endl;
  for (int i=1; i<numLevels; ++i)
  {
    //printf("%d\n", i);
    // Each level is a blurred and downsampled version of the level above it (REDUCE operation)
    //blurredImages[i] = downsample2x(separable_convolution_filter(blurredImages[i-1], kernel,kernel));
    blurredImages[i] = subsample(separable_convolution_filter(blurredImages[i-1], kernel,kernel), 2);
    //std::cout << "blurred: " << blurredImages[i].cols() << " " << blurredImages[i].rows() << std::endl;
    
    std::stringstream path;
    path << "/home/smcmich1/data/demBlendTest/pyrBlur" << i << ".tiff";
    write_image(path.str(), blurredImages[i]);
  }
  
  // Now set up the Laplacian representation
  printf("Building Laplacian pyramid\n");
  for (int i=0; i<numLevels-1; ++i)
  {
    //printf("%d\n", i);
    // Each level is diff between current gaussian level and EXPAND operation of level below it.
    //std::cout << "upsample: " << upsample2x(blurredImages[i+1]).cols() << " " << upsample2x(blurredImages[i+1]).rows() << std::endl;
    //std::cout << "upsample: " << separable_convolution_filter(upsample2x(blurredImages[i+1]), kernel,kernel).cols() << " " << separable_convolution_filter(upsample2x(blurredImages[i+1]), kernel,kernel).rows() << std::endl;
    
    pyramid[i] = blurredImages[i] - 4.0*separable_convolution_filter(upsample2x(blurredImages[i+1]), kernel,kernel);
    
    std::stringstream path;
    path << "/home/smcmich1/data/demBlendTest/pyrStart" << i << ".tiff";
    write_image(path.str(), pyramid[i]);
    
  }
  pyramid[numLevels-1] = blurredImages[numLevels-1];
  
  std::stringstream path2;
  path2 << "/home/smcmich1/data/demBlendTest/pyrStart" << numLevels-1 << ".tiff";
  write_image(path2.str(), pyramid[numLevels-1]);

  // The bottom pyramid level contains the last subsampled image.

  printf("Finished building Laplacian pyramid\n");

  return true;     
}


//TODO: May need to add vw "magic" to make inputs more flexible.
template <typename T>
bool collapseLaplacianImagePyramid(const std::vector<ImageView<T> > &pyramid,
                                      ImageView<T> &outputImage)
{
  // Set up the gaussian smoothing kernel
  // - Copied from stereo PyramidCorrelationView, experiment with this.
  // - This is a 5x5 seperable kernel.
  std::vector<typename DefaultKernelT<T>::type > kernel(5);
  kernel[0] = kernel[4] = 1.0/16.0;
  kernel[1] = kernel[3] = 4.0/16.0;
  kernel[2] = 6.0/16.0;


  // Initialize temporary storage for blurred images
  // TODO: Any way to avoid allocating these buffers?
  int numLevels = (int)pyramid.size();
  std::vector<ImageView<T> > blurredImages(numLevels);

  printf("Starting pyramid collapse\n");
  blurredImages[numLevels-1] = pyramid[numLevels-1]; //TODO: Avoid copies

  std::stringstream path;
  path << "/home/smcmich1/data/demBlendTest/pyrEnd" << numLevels-1 << ".tiff";
  write_image(path.str(), pyramid[numLevels-1]);

  for (int i=numLevels-2; i>=0; --i) 
  {
    //printf("%d\n", i);
    //printf("Reconstructing upwards\n");
    blurredImages[i] = pyramid[i] + 4.0*separable_convolution_filter(upsample2x(blurredImages[i+1]), kernel,kernel);
    std::stringstream path;
    path << "/home/smcmich1/data/demBlendTest/pyrEnd" << i << ".tiff";
    write_image(path.str(), pyramid[i]);
    
  }
  outputImage = copy(blurredImages[0]); // TODO: Avoid copies
  printf("Finished collapsing pyramid\n");
  
  //exit(1);
    
  return true;     
}


/* Test code for Laplacian pyramid
int main(int argc, char *argv[])
{
  
  ImageView<PixelGray<float32> > inputImage, outputImage, diff;
  
  std::vector<ImageView<PixelGray<float32> > > pyramid;
  
  printf("Reading input image\n");
  read_image(inputImage, "/home/smcmich1/data/pyramid/in.tiff");
  write_image("/home/smcmich1/data/pyramid/inBack.tiff", inputImage);
  
  printf("Building pyramid\n");
  buildLaplacianImagePyramid(inputImage, 3, pyramid);
  
  printf("Writing pyramid\n");
  write_image("/home/smcmich1/data/pyramid/level0.tiff", pyramid[0]);
  write_image("/home/smcmich1/data/pyramid/level1.tiff", pyramid[1]);
  write_image("/home/smcmich1/data/pyramid/level2.tiff", pyramid[2]);
  
  
  printf("Collapsing pyramid\n");
  collapseLaplacianImagePyramid(pyramid, outputImage);
  
  printf("Writing output image\n");
  write_image("/home/smcmich1/data/pyramid/out.tiff", outputImage);
  
  diff = inputImage - outputImage;
  write_image("/home/smcmich1/data/pyramid/diff.tiff", diff);
  
  return 0;
}
*/

/// Return an estimate of the georeference resolution
/// - Output format is [dLon, dLat]
/// - Can't use meters since not all GeoReference objects go to meters.
Vector2 getGeorefResolution(const vw::cartography::GeoReference &georef)
{
  Vector2 llo, llx, lly;
  llo = georef.pixel_to_lonlat(Vector2(0,0));
  llx = georef.pixel_to_lonlat(Vector2(1,0));
  lly = georef.pixel_to_lonlat(Vector2(0,1));
  return Vector2(llx[0] - llo[0], lly[1] - llo[1]);
}


class DemMosaicView: public ImageViewBase<DemMosaicView>{
  
  /// Input pixels are value and alpha channel
  typedef           PixelGrayA<float>   DemPixelAlphaT;
  typedef PixelMask<PixelGray <float> > DemPixelMaskT;

  int    m_cols, m_rows;
  bool   m_use_no_weights;
  bool   m_stack_dems;
  float m_out_nodata_value;
  //vector<ModelParams              > const& m_modelParamsArray;
  vector<DiskImageView<DemPixelAlphaT> > const& m_images;
  vector<cartography::GeoReference> const& m_georefs; 
  vector<float>                    const& m_nodata_values;
  cartography::GeoReference m_out_georef;


public:
  DemMosaicView(int cols, int rows, bool use_no_weights, bool stack_dems, float out_nodata_value,
                //vector<ModelParams              > const& modelParamsArray,
                vector<DiskImageView<DemPixelAlphaT> > const& images,
                vector<cartography::GeoReference> const& georefs,
                vector<float                   > const& nodata_values,
                cartography::GeoReference         const& out_georef
                ):
    m_cols(cols), m_rows(rows),
    m_use_no_weights  (use_no_weights), 
    m_stack_dems      (stack_dems),
    m_out_nodata_value(out_nodata_value),
    //m_modelParamsArray(modelParamsArray),
    m_images          (images), 
    m_georefs         (georefs),
    m_nodata_values   (nodata_values),
    m_out_georef      (out_georef)
    {}
  
  typedef float                                  pixel_type;
  typedef pixel_type                             result_type;
  typedef PixelMask<float>                       masked_pixel_type;
  typedef ProceduralPixelAccessor<DemMosaicView> pixel_accessor;
  
  inline int32 cols  () const { return m_cols; }
  inline int32 rows  () const { return m_rows; }
  inline int32 planes() const { return 1; }
  
  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  inline pixel_type operator()( double i, double j, int32 p  = 0 ) const {
    vw_throw(NoImplErr() << "DemMosaicView::operator()(...) is not implemented");
    return pixel_type();
  }

  typedef CropView<ImageView<pixel_type> > prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {   
    
    vw_out() << "------Rasterizing tile of size: " << bbox << "\n";
    // The tile size is always a power of two
    
    
    // Determine the number of pyramid layers
    // - This should be an input parameter!
    // - This is based on the output resolution.
    // - The interpolation calls will take care of input 
    //   DEMs which are not the same resolution.
    const int numPyramidLayers = 3; //TODO: Make a user input!
    
    printf("Setting up accumulators...\n");
    // Set up pyramid accumulators, one for value and one for weight.
    // - Thi size of each pyramid layer will be a power of two.
    size_t layer_width  = bbox.width();
    size_t layer_height = bbox.height();
    std::vector<ImageView<pixel_type> > tile_pyramid  (numPyramidLayers);
    std::vector<ImageView<pixel_type> > weight_pyramid(numPyramidLayers);
    for (int l=0; l<numPyramidLayers; ++l) {
      tile_pyramid  [l].set_size(layer_width, layer_height);
      weight_pyramid[l].set_size(layer_width, layer_height);
      fill(tile_pyramid  [l], m_out_nodata_value);     
      fill(weight_pyramid[l], 0.0);
      layer_width  /= 2;
      layer_height /= 2;
    }
    
    ImageView<pixel_type> output_tile(bbox.width(), bbox.height()); // The output tile!

    // Loop through all input DEMs
    const int numDems = (int)m_images.size();
    for (int dem_iter = 0; dem_iter < numDems; dem_iter++){
      vw_out() << "Starting DEM " << dem_iter << std::endl;


      // Select the georeference and image data for the current DEM
      cartography::GeoReference    georef        = m_georefs[dem_iter];
      ImageViewRef<DemPixelAlphaT> curr_disk_dem = m_images [dem_iter];
      double                       nodata_value  = m_nodata_values[dem_iter];

      
      // Get the tile corner lat/lon locations using the output DEM georeference
      BBox2i bbox_bigger = bbox;
      bbox_bigger.expand(bbox.width()/10); //TODO: 
      Vector2 minLonLat = m_out_georef.pixel_to_lonlat(bbox_bigger.min());
      Vector2 maxLonLat = m_out_georef.pixel_to_lonlat(bbox_bigger.max());
      
      // Determine the pixels in the current DEM that those tile corners correspond to.
      Vector2 minPixel = georef.lonlat_to_pixel(minLonLat);
      Vector2 maxPixel = georef.lonlat_to_pixel(maxLonLat);
      
      // Make sure that the max and min pixel coordinates are not swapped
      if (minPixel[0] > maxPixel[0]) std::swap(minPixel[0], maxPixel[0]);
      if (minPixel[1] > maxPixel[1]) std::swap(minPixel[1], maxPixel[1]);
      minPixel = floor(minPixel); // Round outwards
      maxPixel = ceil (maxPixel);
           
      // Set up the bounding box in the current DEM corresponding the the current output tile
      BBox2i curr_box(minPixel[0], minPixel[1], maxPixel[0] - minPixel[0], maxPixel[1] - minPixel[1]);

      //double dem_width_scale  = bbox.width()  / curr_box.width();
      //double dem_height_scale = bbox.height() / curr_box.height();
      //printf("dem_width_scale  = %lf\n", dem_width_scale);
      //printf("dem_height_scale = %lf\n", dem_height_scale);
      
      //TODO: Do something to make sure outer DEM boundary is sampled!

      curr_box.expand(BilinearInterpolation::pixel_buffer + 1); // Leave an interpolation buffer     
      curr_box.crop(bounding_box(curr_disk_dem));               // Restrict to the DEM size in pixels.
      if (curr_box.empty()) {
        vw_out() << "Skipping DEM, no overlap with tile " << std::endl;
        continue; // If no pixels for this DEM fall in this tile, move to the next tile.      
      }
      vw_out() << "Using BBox " << curr_box << std::endl;
      
      printf("Loading DEM...\n");
      
      // Just load the DEM from disk to interpolate from
      ImageView<DemPixelAlphaT> raw_disk_dem = crop(curr_disk_dem, curr_box);
      // Set up an interpolation wrapper around the loaded DEM image.
      ImageViewRef<DemPixelAlphaT> interp_raw_disk_dem
        = interpolate(raw_disk_dem, BilinearInterpolation(), ConstantEdgeExtension());
      
      printf("Resampling DEM to output tile size...\n");
      // Initialize the resampled DEM to the desired size
      ImageView<DemPixelAlphaT> curr_dem(bbox.width(), bbox.height());
      fill(curr_dem, DemPixelAlphaT(0, 0)); // Initialize to all invalid pixels
      
      //TODO: Write a function to do this!
      //TODO: Should loop be replaced by a functor or something?
      // Loop through each pixel and interpolate from the raw input DEM.
      for (int c = 0; c < bbox.width(); c++){ 
        for (int r = 0; r < bbox.height(); r++){
                                
            // Convert from tile pixel coords to output pixel coords
            Vector2 out_pix(c +  bbox.min().x(), r +  bbox.min().y());
            
            // Get the lon-lat location of this pixel
            Vector2 thisLonLat = m_out_georef.pixel_to_lonlat(out_pix);
            
            // Get the corresponding pixel in the input DEM
            Vector2 in_pix = georef.lonlat_to_pixel(thisLonLat);
            
            // Convert from absolute pixel coords to loaded DEM coords
            double x = in_pix[0] - curr_box.min().x();
            double y = in_pix[1] - curr_box.min().y();
            
            // Logic note: x/y are doubles.
            if ((x < 0) || (x > interp_raw_disk_dem.cols()-1) ||
                (y < 0) || (y > interp_raw_disk_dem.rows()-1) ) 
              continue; // Don't interpolate outside the loaded data, just leave this pixel invalid.
              
            // If we have weights of 0, that means there are invalid pixels, so skip this point.
            int i = (int)floor(x); // Round to a whole pixel because we don't
            int j = (int)floor(y); //  want to interpolate the alpha channel here.
            if (interp_raw_disk_dem(i,   j  ).a() <= 0 ||
                interp_raw_disk_dem(i+1, j  ).a() <= 0 ||
                interp_raw_disk_dem(i,   j+1).a() <= 0 ||
                interp_raw_disk_dem(i+1, j+1).a() <= 0 ) 
              continue;
            
            // Now interpolate the input DEM value and the alpha value and assign to resampled DEM.
            curr_dem(c,r) = interp_raw_disk_dem(x, y);
        } // End row loop
      } // End col loop
      
      // curr_dem is now initialized with the input DEM resampled to the output resolution!
      
      //output_tile = select_channel(curr_dem,0); //DEBUG!!!!
      
/*                                                 
      std::stringstream path;
      path << "/home/smcmich1/data/demBlendTest/resize_dem" << dem_iter << ".tiff";
      write_image(path.str(), curr_dem);  
                                                
      // DEBUG
      ImageView<DemPixelAlphaT> ext_dem = crop(
                                               edge_extend(curr_disk_dem, 
                                                           ValueEdgeExtension<DemPixelAlphaT>(DemPixelAlphaT(0,0))
                                                          ),
                                               curr_box
                                              );
      std::stringstream path2;
      path2 << "/home/smcmich1/data/demBlendTest/ext_dem" << dem_iter << ".tiff";
      write_image(path2.str(), ext_dem);  
                                                 
      //output_tile = select_channel(curr_dem,0);  //DEBUG: Just copy the input DEM to the output tile!
*/      
                        
      // Set up an interpolation wrapper around the loaded DEM image.
      ImageViewRef<DemPixelAlphaT> interp_dem
        = interpolate(curr_dem, BilinearInterpolation(), ConstantEdgeExtension());
        
      // Generate another copy of the input image with the alpha 
      //  channel replaced with a pixel mask.
      ImageView<DemPixelMaskT> curr_dem_masked = alpha_to_mask(curr_dem);
      ImageViewRef<DemPixelMaskT> interp_dem_masked
        = interpolate(curr_dem_masked, BilinearInterpolation(), ConstantEdgeExtension());        
      
      
      printf("Building pyramid...\n");
      // Generate laplacian pyramid for the current DEM 
      // - This replaces the alpha channel with a pixel mask.
      // - TODO: May need to add some special handling of invalid pixels so that they do not
      //         erode away too much of the lower pyramid levels.
      std::vector<ImageView<DemPixelMaskT> > demPyramid;
      if (!buildLaplacianImagePyramid(interp_dem_masked, numPyramidLayers, demPyramid)) {
        vw_throw(InputErr() << "DemMosaicView::Failed to create pyramid for input DEM!");
      }
      
      
      // Loop through all pyramid levels
      // - Layer zero is the highest resolution, down to the lowest at numPyramidLayers-1.
      for (int layer=0; layer<numPyramidLayers; ++layer) {
      
        vw_out() << "Starting layer " << layer << std::endl;
        
        int    layer_width   = tile_pyramid[layer].cols();
        int    layer_height  = tile_pyramid[layer].rows();
        double res_reduction = pow(2,layer);
        
        // Generate alpha scaling
        // - Input alpha map (not a pyramid) has the largest alpha spread, used for the lowest res level.
        // - Each res level above that should use half the blending distance.
        //   This is achieved by multiplying the alpha by a factor of two.
        double alpha_scale  = pow(1.4, numPyramidLayers -1 -layer);

        printf("res_reduction = %lf\n", res_reduction);
        printf("alpha_scale   = %lf\n", alpha_scale);

        // Loop through all pixels in the output pyramid layer
        for (int c = 0; c < layer_width; c++){ 
          for (int r = 0; r < layer_height; r++){

            // Convert to the full-res pixel coordinate on the current tile
            double ct = c * res_reduction;
            double rt = r * res_reduction;
            
            // Now interpolate the input DEM value and the alpha value.
            float val = demPyramid[layer](c, r); // Value from the pyramid level
            float wt  = interp_dem(ct, rt).a();  // The alpha is always computed at full res
            
            // Apply the alpha scaling factor to the weight
            wt *= alpha_scale;
          

            if (wt <= 0.0) // Double check that the alpha value is valid.
              continue;
            if (wt > 1.0) // The weight saturates at 1.0
              wt = 1.0;
            
            //TODO: Retain these options?
            //if (m_use_no_weights || m_stack_dems) 
            //  wt = 1.0; // If either is true the weight must be 0 or 1.
              
            //printf("c, r = %d, %d\n", c, r);
              
            float current_val    = tile_pyramid  [layer](c, r);
            float current_weight = weight_pyramid[layer](c, r);
            if ( current_val == m_out_nodata_value || isnan(current_val) ) {
              tile_pyramid[layer](c, r) = 0; // Clear nodata values
              current_val = 0;
              //printf("Resetting tile value\n");
            }
            
              
            if (dem_iter == (numDems-1)) { // If on the last (priority) DEM
              if (wt >= 1.0){ // Ignore previous values
                tile_pyramid  [layer](c, r) = val;
                weight_pyramid[layer](c, r) = 1.0;
              } else { // Treat other values as weight (1-wt)
                float inv_wt = 1.0 - wt;
                float norm_val = current_val / current_weight; // Unweighted old sum
                tile_pyramid  [layer](c, r) = inv_wt*norm_val + wt*val;
                weight_pyramid[layer](c, r) = 1.0; // Sum of weights is now 1.0
                /*
                printf("c,r = %d, %d\n", c, r);
                printf("Current val = %f\n", current_val);
                printf("Current wt  = %f\n", current_weight);
                printf("Inv weight  = %f\n", inv_wt);
                printf("Norm val   = %f\n", norm_val);
                printf("Weight  = %f\n", wt);
                printf("Dem val = %f\n", val);
                printf("New val = %f\n", tile_pyramid[layer](c, r));
                */
              }
            }else{ // Normal operation
              // Combine the values
              tile_pyramid  [layer](c, r) += wt*val;
              weight_pyramid[layer](c, r) += wt;
            }
            
            //if (!m_stack_dems){ // Combine the values
            //  tile(c, r) += wt*val;
            //  weights(c, r) += wt;
            //}else{              // Use just the last value
            //  tile(c, r) = wt*val;
            //  weights(c, r) = wt;
            //} // End stack option
            //

          } // End loop through rows
        } // End loop through cols
        
      } // End loop through layers

    } // end iterating over DEMs


    printf("Normalizing values according to weights...\n");

    // Loop through all pyramid levels
    // - Layer zero is the highest resolution, down to the lowest at numPyramidLayers-1.
    for (int layer=0; layer<numPyramidLayers; ++layer) {
             
      // Divide by the weights
      int num_valid_pixels = 0;
      for (int c = 0; c < tile_pyramid[layer].cols(); c++){
        for (int r = 0; r < tile_pyramid[layer].rows(); r++){
          float thisWeight = weight_pyramid[layer](c, r);
          if ( thisWeight > 0 ){
            tile_pyramid[layer](c,r) = tile_pyramid[layer](c,r) / thisWeight;
            num_valid_pixels++;
          }
          
        } // End row loop
      } // End col loop
    } // End layer loop


    printf("Collapsing pyramid\n");
    // Collapse the pyramid
    collapseLaplacianImagePyramid(tile_pyramid, output_tile);

//     vw_out() << "Num valid pixels in " << bbox  << ' ' << num_valid_pixels
//              << std::endl;
    return prerasterize_type(output_tile, -bbox.min().x(), -bbox.min().y(),
                             cols(), rows() );
  }

  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

int main( int argc, char *argv[] ) {

  for (int s = 0; s < argc; s++) vw_out() << argv[s] << " ";
  vw_out() << endl;

  try{

    bool use_no_weights = false, stack_dems = false;
    string dem_list_file, out_dem_dir;
    double mpp = 0.0, tr = 0.0;
    float out_nodata_value = numeric_limits<float>::quiet_NaN();
    int num_threads = 0, tile_size = -1, tile_index = -1;
    po::options_description general_options("Options");
    general_options.add_options()
      ("dem-list-file,l", po::value<string>(&dem_list_file),
       "List of DEM files to mosaic, one per line.")
      ("tr", po::value(&tr)->default_value(0.0),
       "Output DEM resolution in target georeferenced units per pixel. If not specified, use the same resolution as the first DEM to be mosaicked.")
      ("mpp", po::value<double>(&mpp)->default_value(0.0),
       "Output DEM resolution in meters per pixel (at the equator). If not specified, use the same resolution as the first DEM to be mosaicked.")
      ("tile-size", po::value<int>(&tile_size)->default_value(1000000),
       "The maximum size of output DEM tile files to write, in pixels.")
      ("tile-index", po::value<int>(&tile_index)->default_value(0),
       "The index of the tile to save in the list of tiles (starting from zero). If called with given tile size, and no tile index, the tool will print out how many tiles are there.")
      ("output-dem-dir,o", po::value<string>(&out_dem_dir),
       "The directory in which to save the DEM tiles.")
      ("output-nodata-value", po::value<float>(&out_nodata_value),
       "No-data value to use on output.")
      ("use-no-weights", po::bool_switch(&use_no_weights)->default_value(false),
       "Average the DEMs to mosaic without using weights (the result is not as smooth).")
      ("stack-dems", po::bool_switch(&stack_dems)->default_value(false),
       "Stack each DEM on top of the previous one, without attempting to combine their values (the result is less smooth).")
      ("threads", po::value<int>(&num_threads),
       "Number of threads to use.")
      ("help,h", "Display this help message.");
    ;

    // Parse options
    po::options_description hidden_options("");
    po::options_description options("Allowed Options");
    options.add(general_options).add(hidden_options);
    po::positional_options_description p;
    ostringstream usage;
    usage << "A tool for mosaicking DEMs without seams.\n\n";
    usage << general_options << endl;
    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).options(options)
               .positional(p).run(), vm );
    po::notify( vm );
    if  (vm.count("help") || argc <= 1){
      vw_out() << usage.str() << "\n";
      return 1;
    }

    // Error checking
    if (dem_list_file == "")
      vw_throw(ArgumentErr() << "No list of DEMs was specified.\n");
    if (mpp > 0.0 && tr > 0.0)
      vw_throw(ArgumentErr() << "Just one of the --mpp and --tr options needs to be set.\n");
    if (out_dem_dir == "")
      vw_throw(ArgumentErr() << "No output DEM directory was specified.\n");
    if (num_threads == 0)
      vw_throw(ArgumentErr() << "The number of threads must be set and "
               << "positive.\n");
    if (tile_size <= 0)
      vw_throw(ArgumentErr() << "The size of a tile in pixels must "
               << "be set and positive.\n");
    if (tile_index < 0)
      vw_throw(ArgumentErr() << "The index of the tile to save must be set "
               << "and non-negative.\n");
    if (use_no_weights && stack_dems){
      vw_throw(ArgumentErr() << "At most one of --use-no-weights and --stack-dems must be set.\n");
    }
    
    // Read the DEMs to mosaic
    vector<string> dem_files;
    ifstream is(dem_list_file.c_str());
    string file;
    while (is >> file) dem_files.push_back(file);
    if (dem_files.empty())
      vw_throw(ArgumentErr() << "No DEM files to mosaic.\n");
    is.close();
    
    // Read nodata from first DEM, unless the user chooses to specify it.
    if (!vm.count("output-nodata-value")){
      DiskImageResourceGDAL in_rsrc(dem_files[0]);
      if ( in_rsrc.has_nodata_read() ) out_nodata_value = in_rsrc.nodata_read();
    }
    
    // Find the lon-lat bounding box of all DEMs
    double big = numeric_limits<double>::max();
    Vector4 ll_bbox(big, -big, big, -big);
    for (int dem_iter = 0; dem_iter < (int)dem_files.size(); dem_iter++){
      std::string curr_file = dem_files[dem_iter];
      Vector4 corners = getImageCorners(curr_file);
      ll_bbox[0] = std::min(ll_bbox[0], corners[0]);
      ll_bbox[1] = std::max(ll_bbox[1], corners[1]);
      ll_bbox[2] = std::min(ll_bbox[2], corners[2]);
      ll_bbox[3] = std::max(ll_bbox[3], corners[3]);
    }

    // Form the georef. The georef of the first DEM is used as initial guess.
    cartography::GeoReference out_georef;
    bool is_good = read_georeference(out_georef, dem_files[0]);
    if (!is_good){
      std::cerr << "No georeference found in " << dem_files[0] << std::endl;
      exit(1);
    }
    double spacing = tr;
    if (mpp > 0.0){
      // First convert meters per pixel to degrees per pixel using
      // the datum radius.
      spacing = 360.0*mpp/( 2*M_PI*out_georef.datum().semi_major_axis() );
      // Next, wipe the georef altogether, as we will use lon-lat.
      out_georef.set_geographic();
    }

    // Use desired spacing if user-specified
    if (spacing > 0.0){
      Matrix<double,3,3> transform = out_georef.transform();
      transform.set_identity();
      transform(0, 0) = spacing;
      transform(1, 1) = -spacing;
      out_georef.set_transform(transform);
    }
    
    // Set the lower-left corner
    Vector2 beg_pix = out_georef.lonlat_to_pixel(Vector2(ll_bbox[0], ll_bbox[3]));
    out_georef = crop(out_georef, beg_pix[0], beg_pix[1]); // Reset georef to start at this pixel
    std::cout << "Output georeference:\n" << out_georef << std::endl;
    
    // Image size
    Vector2 end_pix = out_georef.lonlat_to_pixel(Vector2(ll_bbox[1], ll_bbox[2]));
    int cols = (int)ceil(end_pix[0]); // Georef starts at zero and goes to these pixels
    int rows = (int)ceil(end_pix[1]);

    //TODO: Force tile_size to be multiple of internal tile size!
    //TODO: Pad to multiple of external tile size!
    int pad_cols = cols % 256;
    int pad_rows = cols % 256;
    if (pad_cols > 0)
      cols += (256-pad_cols);
    if (pad_rows > 0)
      rows += (256-pad_rows);

    typedef PixelGrayA<float> DemPixelAlphaT; //TODO: Define this somewhere else?

    // Compute the weights, and store the no-data values, pointers
    // to images, and georeferences (for speed).
    //vector<ModelParams> modelParamsArray;
    vector<float                        > nodata_values;
    vector<DiskImageView<DemPixelAlphaT> > images;
    vector<cartography::GeoReference     > georefs;
    for (int dem_iter = 0; dem_iter < (int)dem_files.size(); dem_iter++){ // Loop through all DEMs
      //modelParamsArray.push_back(ModelParams());
      //modelParamsArray[dem_iter].inputFilename = dem_files[dem_iter];
      images.push_back(DiskImageView<DemPixelAlphaT>( dem_files[dem_iter] ));

      cartography::GeoReference geo;
      bool is_good = read_georeference(geo, dem_files[dem_iter]);
      if (!is_good){
        std::cerr << "No georeference found in " << dem_files[dem_iter]
                  << std::endl;
        exit(1);
      }
      georefs.push_back(geo);
      
      float curr_nodata_value = out_nodata_value;
      DiskImageResourceGDAL in_rsrc(dem_files[dem_iter]);
      if ( in_rsrc.has_nodata_read() ) curr_nodata_value = in_rsrc.nodata_read();
      nodata_values.push_back(curr_nodata_value);

      boost::scoped_ptr<vw::SrcImageResource> src(vw::DiskImageResource::open(dem_files[dem_iter]));
      int num_channels = src->channels();
      int num_planes   = src->planes();
      if (num_channels*num_planes != 2){
        std::cerr << "Need to have an alpha channel with the grassfire weights in "
                  << dem_files[dem_iter] << std::endl;
        exit(1);
      }

    }


    // Form the mosaic and write it to disk
    vw_out()<< "The size of the mosaic is " << cols << " x " << rows
            << " pixels.\n";

    int num_tiles_x = (int)ceil((double)cols/double(tile_size));
    if (num_tiles_x <= 0) num_tiles_x = 1;
    int num_tiles_y = (int)ceil((double)rows/double(tile_size));
    if (num_tiles_y <= 0) num_tiles_y = 1;
    int num_tiles = num_tiles_x*num_tiles_y;
    vw_out() << "Number of tiles: " << num_tiles_x << " x "
             << num_tiles_y << " = " << num_tiles << std::endl;

    if (tile_index < 0 || tile_index >= num_tiles){
      vw_out() << "Tile with index: " << tile_index << " is out of bounds."
               << std::endl;
      return 0;
    }
    int tile_index_y = tile_index / num_tiles_y;
    int tile_index_x = tile_index - tile_index_y*num_tiles_y;
    BBox2i tile_box(tile_index_x*tile_size, tile_index_y*tile_size,
                    tile_size, tile_size);
    tile_box.crop(BBox2i(0, 0, cols, rows));
    ostringstream os; os << out_dem_dir << "/tile-" << tile_index << ".tif";
    std::string dem_tile = os.str();
    create_out_dir(dem_tile);
    vw::vw_settings().set_default_num_threads(num_threads);
    DiskImageResourceGDAL::Options gdal_options;
    gdal_options["COMPRESS"] = "LZW";
    ImageViewRef<float> out_dem
      = crop(DemMosaicView(cols, rows, use_no_weights, stack_dems, out_nodata_value, //modelParamsArray, 
                           images, georefs, nodata_values, out_georef),
             tile_box);            
    vw_out() << "Writing: " << dem_tile << std::endl;
    DiskImageResourceGDAL rsrc(dem_tile, out_dem.format(), Vector2i(256, 256),
                               gdal_options);
    if (!isnan(out_nodata_value))
      rsrc.set_nodata_write(out_nodata_value);
    cartography::GeoReference crop_georef
      = crop(out_georef, tile_box.min().x(), tile_box.min().y());
    
    write_georeference(rsrc, crop_georef);
    block_write_image(rsrc, out_dem, TerminalProgressCallback("{Core}",
                                                              "Processing:"));
    
  } catch ( const exception& e ) {
    cerr << "\n\nError: " << e.what() <<  endl;
    return 1;                                     
  }

  return 0;
}


