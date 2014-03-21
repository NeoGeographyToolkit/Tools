NGT Tools
=====

This repository contains small programs which may be useful while working with Stereo Pipeline.  They are typically not as fully tested as the tools packaged with Stero Pipeline but could be cleaned up and incorporated into Stereo Pipeline if there is enough interest.  In order to use these tools you must first have Vision Workbench and Stereo Pipeline installed using BinaryBuilder.



--- Detailed build instructions ---

Build the tools using CMake (version >= 2.8) using the following series of commands:

git clone https://github.com/NeoGeographyToolkit/Tools.git
cd Tools
mkdir build
cd build
cmake -D BASESYSTEM_INSTALL_DIR=<path to your base sytem install directory> -D VISIONWORKBENCH_INSTALL_DIR=<path to your vision workbench install directory> -D STEREOPIPELINE_INSTALL_DIR=<path to your stereo pipeline install directory> -D CMAKE_BUILD_TYPE=<probably Release>..




