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

#include <vw/Core.h>
#include <vw/FileIO.h>
#include <vw/Image.h>
#include <vw/Cartography.h>


#include <asp/Core/Macros.h>
#include <asp/Core/Common.h>
#include <asp/Sessions/CameraModelLoader.h>
#include <asp/asp_config.h>
#include <asp/Core/BundleAdjustUtils.h>

using namespace vw;
using namespace asp;
/*
Vector2 demPixToCamPix(Vector2i const& dem_pixel,
                      boost::shared_ptr<camera::CameraModel> const& camera_model,
                      ImageViewRef<PMaskT> const& dem,
                      GeoReference const &dem_georef)
{
  Vector2 lonlat = dem_georef.point_to_lonlat(dem_georef.pixel_to_point(dem_pixel));
  //vw_out() << "lonlat = " << lonlat << std::endl;
  PMaskT height = dem(dem_pixel[0], dem_pixel[1]);
  Vector3 xyz = dem_georef.datum().geodetic_to_cartesian
                    (Vector3(lonlat[0], lonlat[1], height.child()));
  //vw_out() << "xyz = " << xyz << std::endl;
  // Throws if the projection fails ???
  Vector2i camera_pixel = camera_model->point_to_pixel(xyz);
  //vw_out() << "camera_pixel = " << camera_pixel << std::endl;
  return camera_pixel;
}
*/
/// Print out some percentiles from one small vector
void printStatistics(const std::vector<double> &vec)
{
  std::vector<double> sorted(vec);
  std::sort(sorted.begin(), sorted.end());
  
  std::vector<double> percentilesToCompute;
  for (int i=1; i<20; ++i)
    percentilesToCompute.push_back((double)i / 20.0);

  for (size_t i=0; i<percentilesToCompute.size(); ++i)
  { 
    double percentile = percentilesToCompute[i];
    int    index      = percentile * vec.size();
    double value      = sorted[index];
    printf("Percentile %lf = %lf\n", percentile, value);
  }
}

typedef boost::shared_ptr<vw::camera::CameraModel> CameraModelPtr;
typedef boost::shared_ptr<asp::RPCModel> RPCModelPtr;


RPCModelPtr loadRpcModel(const std::string &path, vw::BBox3 &validBounds, vw::Vector2i &imageSize)
{
  RPCXML rpc_xml; // This is for reading XML files
  rpc_xml.read_from_file(path);
  RPCModel* rpc_model = new RPCModel(*rpc_xml.rpc_ptr()); // Copy the value

  // We also need the valid bounding box
  validBounds = rpc_xml.get_lon_lat_height_box();

  imageSize = xml_image_size(path);

  std::cout << "valid bounds = " << validBounds << std::endl;
  std::cout << "Image size = " << imageSize << std::endl;

  return boost::shared_ptr<asp::RPCModel>(rpc_model);
}



int main( int argc, char* argv[] ) {

  if (argc < 5)
  {
    printf("Usage: sensorModelTest <dg model path> <rpc model path 1> <rpc model path 2> <output prefix>\n");
    return -1;
  }
    
  std::string dgCameraFile = argv[1];
  std::string theirRpcFile = argv[2];
  std::string ourRpcFile   = argv[3];
  std::string outputPrefix = argv[4];

  std::ofstream theirErrorFile(std::string(outputPrefix + "_their_rpc_error.csv"  ).c_str());
  std::ofstream ourErrorFile  (std::string(outputPrefix + "_our_rpc_error.csv"    ).c_str());
  std::ofstream rpcErrorFile  (std::string(outputPrefix + "_rpc_rpc_error.csv"    ).c_str());

  std::ofstream allOutputFile (std::string(outputPrefix + "_all_model_outputs.csv").c_str());

  // Load the DG sensor model
  printf("\nLoading DG sensor model...\n");
  asp::CameraModelLoader cameraModelLoader;
  CameraModelPtr dgCameraModel = cameraModelLoader.load_dg_camera_model(dgCameraFile);

  // Load the existing RPC sensor model
  printf("\nLoading DG's RPC sensor model and valid boundaries...\n");
  vw::BBox3    validBounds;
  vw::Vector2i imageSize;
  RPCModelPtr theirRpcCameraModel = loadRpcModel(theirRpcFile, validBounds, imageSize);
  if (!theirRpcCameraModel)
    return -1;


  vw::cartography::Datum datum = theirRpcCameraModel->datum();

  // Call our solver to compute an estimated RPC model
  // Or just load our computed RPC model.
  printf("\nLoading our solved RPC sensor model...\n");
  RPCModelPtr ourRpcCameraModel = cameraModelLoader.load_rpc_camera_model(ourRpcFile);

  std::vector<double> theirRpcDiffList, ourRpcDiffList, rpcDiffList;

  const bool TEST_RAY_GENERATION = false;


  // Extract the area of interest (AOI) bounds
  double xMin, xMax, yMin, yMax, zMin, zMax; 

  if (!TEST_RAY_GENERATION) // Test point projection
  {

	  xMin = validBounds.min().x();  // These are lon/lat coordinates
	  xMax = validBounds.max().x();
	  yMin = validBounds.min().y();
	  yMax = validBounds.max().y();
    zMin = validBounds.min().z(); // There may only be a single z value among the four points
	  zMax = validBounds.max().z();


    int numPointsX = 100;
    int numPointsY = 100;
    int numPointsZ = 1;
    int totalNumPoints = numPointsX * numPointsY * numPointsZ;
    printf("Will compute error at %d points.\n", totalNumPoints);
    theirRpcDiffList.resize(totalNumPoints);
    ourRpcDiffList.resize(totalNumPoints);
    rpcDiffList.resize(totalNumPoints);

    const std::string errorHeader = "# x, y, z, error\n";
    theirErrorFile << errorHeader;
    ourErrorFile   << errorHeader;
    rpcErrorFile   << errorHeader;

    //allOutputFile << "# For each location: [dg] col, row, [their RPC], [our RPC]\n";
    allOutputFile << numPointsX << ", " << numPointsY << std::endl;

    double xStep = (xMax - xMin) / (double)numPointsX;
    double yStep = (yMax - yMin) / (double)numPointsY;
    //double zStep = (zMax - zMin) / (double)numPointsZ;


    // Evaluate the three models at a range of points throughout the AOI
    printf("Computing projection error...\n");
    int count = 0;
    double x, y, z;
    for (int zi=0; zi<numPointsZ; ++zi){
      for (int yi=0; yi<numPointsY; ++yi){
        for (int xi=0; xi<numPointsX; ++xi){
        

          x = xMin + xStep*xi;
          y = yMin + yStep*yi;
          //y = (yMin + yMax)/2.0;
          z = (zMin + zMax)/2.0;//zMin + zStep*zi;

          //try {

            // Convert from GDC to GCC
            Vector3 gdcCoord(x, y, z); // The input bounds are in lon/lat
            Vector3 gccCoord = datum.geodetic_to_cartesian(gdcCoord);

            // Project to a pixel in all three cameras
            Vector2 pixelDg        = dgCameraModel->point_to_pixel(gccCoord);
            Vector2 pixelRpcTheirs = theirRpcCameraModel->point_to_pixel(gccCoord);
            Vector2 pixelRpcOurs   = ourRpcCameraModel->point_to_pixel(gccCoord);
            /*
            std::cout << std::setprecision(5);
            std::cout << std::endl<< "GDC = "       << gdcCoord << std::endl;
            //std::cout << "GCC = "       << gccCoord << std::endl;
            std::cout << "DG        = " << pixelDg  << std::endl;
            std::cout << "Their RPC = " << pixelRpcTheirs << std::endl;
            //std::cout << "Our RPC   = " << pixelRpcOurs   << std::endl;
  */
            double theirRpcDiff = norm_2(pixelRpcTheirs - pixelDg);
            double ourRpcDiff   = norm_2(pixelRpcOurs   - pixelDg);
            double rpcDiff      = norm_2(pixelRpcOurs   - pixelRpcTheirs);

            /*
            // DEBUG -> Use normalized values
            Vector2 normPixelDg        = vw::math::elem_quot(pixelDg - theirRpcCameraModel->xy_offset(), theirRpcCameraModel->xy_scale());
            Vector2 normPixelRpcTheirs = vw::math::elem_quot(pixelRpcTheirs - theirRpcCameraModel->xy_offset(), theirRpcCameraModel->xy_scale());
            Vector2 normPixelRpcOurs   = vw::math::elem_quot(pixelRpcOurs - theirRpcCameraModel->xy_offset(), theirRpcCameraModel->xy_scale());
            double theirRpcDiff = norm_2(normPixelRpcTheirs - normPixelDg);
            double ourRpcDiff   = norm_2(normPixelRpcOurs   - normPixelDg);
            double rpcDiff      = norm_2(normPixelRpcOurs   - normPixelRpcTheirs);
            */

            //std::cout << "Their diff: " << theirRpcDiff << std::endl;

            //theirErrorFile << theirRpcDiff << std::endl;
            theirErrorFile << x << ", " << y << ", " << z << ", "<< theirRpcDiff << std::endl;
            ourErrorFile   << x << ", " << y << ", " << z << ", "<< ourRpcDiff   << std::endl;
            rpcErrorFile   << x << ", " << y << ", " << z << ", "<< rpcDiff      << std::endl;

            allOutputFile << pixelDg.x()        << ", " << pixelDg.y()        << ", "
                          << pixelRpcTheirs.x() << ", " << pixelRpcTheirs.y() << ", "
                          << pixelRpcOurs.x()   << ", " << pixelRpcOurs.y()   << std::endl;


            theirRpcDiffList[count] = theirRpcDiff;
            ourRpcDiffList  [count] = ourRpcDiff;
            rpcDiffList     [count] = rpcDiff;/*
          } catch (...) // TODO: Should not be any!
          {
            printf("Caught exception!\n");
            theirRpcDiffList[count] = 0;
            ourRpcDiffList  [count] = 0;
            rpcDiffList     [count] = 0;
          }*/

          count++;
        }//x loop
      }//y loop
    }//z loop

  }
  else // Test ray generation
  {
    xMin = 0;
    yMin = 0;
    xMax = imageSize[0] - 1;
    yMax = imageSize[1] - 1;  

    int numPointsX = 100;
    int numPointsY = 100;
    int totalNumPoints = numPointsX * numPointsY;
    printf("Will compute error at %d points.\n", totalNumPoints);
    theirRpcDiffList.resize(totalNumPoints);
    ourRpcDiffList.resize(totalNumPoints);
    rpcDiffList.resize(totalNumPoints);

    const std::string errorHeader = "# col, row, error\n";
    theirErrorFile << errorHeader;
    ourErrorFile   << errorHeader;
    rpcErrorFile   << errorHeader;

    //allOutputFile << "# For each location: [dg] col, row, [their RPC], [our RPC]\n";
    allOutputFile << numPointsX << ", " << numPointsY << std::endl;

    double xStep = (xMax - xMin) / numPointsX;
    double yStep = (yMax - yMin) / numPointsY;


    // Evaluate the three models at a range of points throughout the AOI
    printf("Computing ray generation error...\n");
    int count = 0;
    double x, y;
    for (int yi=0; yi<numPointsY; ++yi){
      for (int xi=0; xi<numPointsX; ++xi){
      
        x = xMin + xStep*xi;
        y = yMin + yStep*yi;
 
        Vector2 pixel(x, y);

        // Generate a ray in all three cameras
        Vector3 rayDg        = dgCameraModel->pixel_to_vector(pixel);
        Vector3 rayRpcTheirs = theirRpcCameraModel->pixel_to_vector(pixel);
        Vector3 rayRpcOurs   = ourRpcCameraModel->pixel_to_vector(pixel);
       
        double theirRpcDiff = norm_2(rayRpcTheirs - rayDg);
        double ourRpcDiff   = norm_2(rayRpcOurs   - rayDg);
        double rpcDiff      = norm_2(rayRpcOurs   - rayRpcTheirs);

        //theirErrorFile << theirRpcDiff << std::endl;
        theirErrorFile << x << ", " << y << ", " << theirRpcDiff << std::endl;
        ourErrorFile   << x << ", " << y << ", " << ourRpcDiff   << std::endl;
        rpcErrorFile   << x << ", " << y << ", " << rpcDiff      << std::endl;

        allOutputFile << rayDg.x()        << ", " << rayDg.y()        << ", " << rayDg.z()        << ", "
                      << rayRpcTheirs.x() << ", " << rayRpcTheirs.y() << ", " << rayRpcTheirs.z() << ", "
                      << rayRpcOurs.x()   << ", " << rayRpcOurs.y()   << ", " << rayRpcOurs.z()   << std::endl;


        theirRpcDiffList[count] = theirRpcDiff;
        ourRpcDiffList  [count] = ourRpcDiff;
        rpcDiffList     [count] = rpcDiff;

        count++;
      }//x loop
    }//y loop


  }



  // Close log files
  theirErrorFile.close();
  ourErrorFile.close();
  rpcErrorFile.close();

  allOutputFile.close();

  // Compute statistics
  printf("Computing statistics:\n");

  printf("\nOur RPC vs DG difference:\n");
  printStatistics(ourRpcDiffList);

  // If we only loaded one RPC model then there is only one set of differences to display
  if (ourRpcFile != theirRpcFile)
  {
    printf("\nDG's RPC vs DG difference:\n");
    printStatistics(theirRpcDiffList);

    printf("\nRPC vs RPC difference:\n");
    printStatistics(rpcDiffList);
  }

  printf("\nFinished!\n");

  return 0;
}









  //std::string testCameraFile = "/home/smcmich1/data/rpc_accuracy/09OCT11191503-P1BS_R1C1-052783426010_01_P001.XML";  
  //std::string ourRpcFile     = "/home/smcmich1/data/rpc_accuracy/09OCT11191503-P1BS_R1C1-052783426010_01_P001.our_rpc2.xml";





/*
	// --> Harcoded from the file 09OCT11191503-P1BS_R1C1-052783426010_01_P001.XML
	xMin = -121.23968436;
	xMax = -121.02862451;
	yMin =   36.93682841;
	yMax =   37.06957892;
  zMin =   44.27;
	zMax =  588.84;
*/
/*
	// --> Harcoded from the file WV01_20130513_102001002290EC00_102001002292A900 (TODO: Rest of them!)
  testCameraFile = "/home/oalexan1/projects/data/stereo_antarctic/WV01_20130513_102001002290EC00_102001002292A900/WV01_13MAY131553060-P1BS-102001002290EC00.xml";
  ourRpcFile = "/home/smcmich1/data/rpc_accuracy/our_rpc_file.xml";
	xMin = -53.90592527;
	xMax = -53.35065155;
	yMin = 72.86664777999999;
	yMax = 73.01384518;
  zMin =  836.1700000000000;
	zMax = 1142.710000000000;
*/
	// --> Harcoded from the file WV02_20130829_10300100277B9A00_1030010026789F00 (TODO: Rest of them!)
  // TODO: when testing our batch output, there is no "their RPC" file, just compare the two results in our file!
  //testCameraFile = "/home/oalexan1/projects/data/stereo_antarctic/WV02_20130829_10300100277B9A00_1030010026789F00/WV02_13AUG291630386-P1BS-10300100277B9A00.xml";
  //ourRpcFile = "/home/smcmich1/data/rpc_accuracy/our_rpc_file.xml";
/*
	xMin = -113.0223235500000;
	xMax = -112.4307609000000;
	yMin = -74.11042309000000;
	yMax = -73.97897283000000;
  zMin = -33.00000000000000 - 100;
	zMax = -33.00000000000000;
*/





