#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __BEGIN_LICENSE__
#  Copyright (c) 2009-2013, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NGT platform is licensed under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance with the
#  License. You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# __END_LICENSE__

"""IrgGeoFunctions.py - Functions for working with different geo-data formats"""

import sys, os, glob, re, shutil, subprocess, string, time, errno

import IrgIsisFunctions


# This function is wrapped here for convenience
# - To put it in one place would require another python functions file.
def getImageSize(imagePath):
    return IrgIsisFunctions.getImageSize(imagePath)



def getGeoTiffBoundingBox(geoTiffPath):
    """Returns (minLon, maxLon, minLat, maxLat) for a geotiff image"""
    
    # Call command line tool silently
    cmd = ['getLatLonBounds', geoTiffPath]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    textOutput, err = p.communicate()

    # Check that the call did not fail
    if (textOutput.find('Failed') >= 0):
        raise Exception('Error: getGeoTiffBoundingBox failed on input image: ' + geoTiffPath)
    
    # Parse the output
    lines = textOutput.split('\n')
    minLat = float( lines[0][ lines[0].find('=')+1 :] )
    maxLat = float( lines[1][ lines[1].find('=')+1 :] )
    minLon = float( lines[2][ lines[2].find('=')+1 :] )
    maxLon = float( lines[3][ lines[3].find('=')+1 :] )
    
    return (minLon, maxLon, minLat, maxLat)


def getImageBoundingBox(filePath):
    """Returns (minLon, maxLon, minLat, maxLat) for a georeferenced image file"""

    extension = os.path.splitext(filePath)[1]
    if '.tif' in extension:
        return getGeoTiffBoundingBox(filePath)
    else:
        return IrgIsisFunctions.getIsisBoundingBox(filePath)
    
    # Any other file types will end up raising some sort of exception

