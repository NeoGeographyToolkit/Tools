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

import IrgIsisFunctions, IrgStringFunctions


# This function is wrapped here for convenience
# - To put it in one place would require another python functions file.
def getImageSize(imagePath):
    return IrgIsisFunctions.getImageSize(imagePath)


def getGdalInfoTagValue(text, tag):
    """Gets the value of a gdal parameter in a [""] tag or None if it is absent."""

    try:
        lineAfterTag = IrgStringFunctions.getLineAfterText(text, tag)
        
        # The remaining line should look like this: ",25],
        commaPos   = lineAfterTag.find(',')
        bracketPos = lineAfterTag.find(']')
        # The value is always returned as a string
        return IrgStringFunctions.convertToFloatIfNumber(lineAfterTag[commaPos+1:bracketPos])
    
    except Exception: # Requested tag was not found
        return None

def getImageGeoInfo(imagePath):
    """Obtains some image geo information from gdalinfo in dictionary format"""
    
    # Call command line tool silently
    cmd = ['gdalinfo', imagePath, '-stats']
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    textOutput, err = p.communicate()
    
    originLine    = IrgStringFunctions.getLineAfterText(textOutput, 'Origin = ')
    pixelSizeLine = IrgStringFunctions.getLineAfterText(textOutput, 'Pixel Size = ')
    
    originVals    = IrgStringFunctions.getNumbersInParentheses(originLine)
    pixelSizeVals = IrgStringFunctions.getNumbersInParentheses(pixelSizeLine)
      
    outputDict = {}
    outputDict['origin']     = originVals
    outputDict['pixel size'] = pixelSizeVals
    
    outputDict['standard_parallel_1'] = getGdalInfoTagValue(textOutput, 'standard_parallel_1')
    outputDict['central meridian']    = getGdalInfoTagValue(textOutput, 'central_meridian')
    
    return outputDict

def getImageStats(imagePath):
    """Obtains some image statistics from gdalinfo"""
    
    if not os.path.exists(imagePath):
        raise Exception('Image file ' + imagePath + ' not found!')
    
    # Call command line tool silently
    cmd = ['gdalinfo', imagePath, '-stats']
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    textOutput, err = p.communicate()
    
    # Statistics are computed seperately for each band
    bandStats = []
    band = 0
    while (True): # Loop until we run out of bands
        # Look for the stats line for this band
        bandString = 'Band ' + str(band+1) + ' Block='
        bandLoc = textOutput.find(bandString)
        if bandLoc < 0:
            return bandStats # Quit if we did not find it
            
        # Now parse out the statistics for this band
        bandMaxStart  = textOutput.find('STATISTICS_MAXIMUM=', bandLoc)
        bandMeanStart = textOutput.find('STATISTICS_MEAN=',    bandLoc)
        bandMinStart  = textOutput.find('STATISTICS_MINIMUM=', bandLoc)
        bandStdStart  = textOutput.find('STATISTICS_STDDEV=',  bandLoc)
               
        bandMax  = IrgStringFunctions.getNumberAfterEqualSign(textOutput, bandMaxStart)
        bandMean = IrgStringFunctions.getNumberAfterEqualSign(textOutput, bandMeanStart)
        bandMin  = IrgStringFunctions.getNumberAfterEqualSign(textOutput, bandMinStart)
        bandStd  = IrgStringFunctions.getNumberAfterEqualSign(textOutput, bandStdStart)
            
        # Add results to the output list
        bandStats.append( (bandMin, bandMax, bandMean, bandStd) )
            
        band = band + 1 # Move to the next band
    


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

