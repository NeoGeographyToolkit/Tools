===========================
===      NGT Tools      ===
===========================

This repository contains small programs which may be useful while working with Stereo Pipeline.  They are typically not as fully tested as the tools packaged with Stero Pipeline but could be cleaned up and incorporated into Stereo Pipeline if there is enough interest.  In order to use these tools you must first have Vision Workbench and Stereo Pipeline installed using BinaryBuilder.



------- Detailed build instructions -------

Build the tools using CMake (version >= 2.8) using the following series of commands:

git clone https://github.com/NeoGeographyToolkit/Tools.git
cd Tools
mkdir build
cd build
cmake -D BASESYSTEM_INSTALL_DIR=<path to your base sytem install directory> -D VISIONWORKBENCH_INSTALL_DIR=<path to your vision workbench install directory> -D STEREOPIPELINE_INSTALL_DIR=<path to your stereo pipeline install directory> -D CMAKE_BUILD_TYPE=<probably Release> ..


------- Summary of files -------

--- Cmake files ---

FindStereoPipeline.csv  = Script for finding a Stereo Pipeline installation.
FindVisionWorkbench.csv = Script for finding a Vision Workbench installation.

--- Python files ---

cub2kml.py = Simple tool to generate a kml file showing the bounding box for an ISIS cube.  
           = Only tested on the moon so far.
           = Could be expanded to work on the Earth and other image types.

IrgAspFunctions.py  = Collection of functions for working with ASP data.
IrgFileFunctions.py = Collection of generic file related functions.
IrgGeoFunctions.py  = Collection of functions for working with geo images.
IrgIsisFunctions.py = Collection of function for working with ISIS data/tools.

--- C++ Files ---

stereo.h = Copy of a file from ASP with many dependencies removed.

dem_mean_std.cc           = Find the mean and std dev DEM of the given input sets of DEMs.
                          = Also output count.tif, showing for each pixel how many times it was encountered in the stack of DEMs.
                
getLatLonBounds.cc        = Print the lat/lon bounding box of an input GeoTif image.

imagestats.cc             = Compute a set of statistics of an input image.

lola_compare.cc           = Compare a set of LOLA points to a GeoTif DEM and compute distance statistics.

maskFromIntersectError.cc = Generate an 8 bit grayscale image showing an ASP stereo image divided into regions by intersection error.

matchBinaryToCsv.cc       = Convert the binary output file of ASP interest point matching functions into a human readable CSV file.

pixelPairsFromStereo.cc   = Generate a CSV file containing matching pixel pairs from a dense ASP stereo output file.

point_to_dem_dist.cc      = Compute the distance from a set of points to a DEM below them.

stereoIpFind.cc           = Search for a set of matching pixel pairs in two images and write them to a binary match file.







