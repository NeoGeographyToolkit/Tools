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

# TODO: This can take a long time due to the stats call!
def getImageGeoInfo(imagePath):
    """Obtains some image geo information from gdalinfo in dictionary format"""
    
    outputDict = {}
    
    # Call command line tool silently
    cmd = ['gdalinfo', imagePath, '-stats']
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    textOutput, err = p.communicate()
    
    # Get the size in pixels
    imageSizeLine = IrgStringFunctions.getLineAfterText(textOutput, 'Size is ')
    sizeVals      = imageSizeLine.split(',')
    outputDict['image_size'] = (int(sizeVals[0]), int(sizeVals[1]))

    # Get origin location and pixel size    
    originLine    = IrgStringFunctions.getLineAfterText(textOutput, 'Origin = ')
    pixelSizeLine = IrgStringFunctions.getLineAfterText(textOutput, 'Pixel Size = ')    
    originVals    = IrgStringFunctions.getNumbersInParentheses(originLine)
    pixelSizeVals = IrgStringFunctions.getNumbersInParentheses(pixelSizeLine)
    outputDict['origin']     = originVals
    outputDict['pixel size'] = pixelSizeVals

    # Get these two values!
    outputDict['standard_parallel_1'] = getGdalInfoTagValue(textOutput, 'standard_parallel_1')
    outputDict['central_meridian']    = getGdalInfoTagValue(textOutput, 'central_meridian')

    
    # Extract this variable which ASP inserts into its point cloud files
    try:
        pointOffsetLine = IrgStringFunctions.getLineAfterText(textOutput, 'POINT_OFFSET=') # Tag name must be synced with C++ code
        offsetValues    = pointOffsetLine.split(' ')
        outputDict['point_offset'] =  (float(offsetValues[0]), float(offsetValues[1]), float(offsetValues[2]))        
    except:
        pass # In most cases this line will not be present

    
    # List of dictionaries per band
    outputDict['band_info'] = []
    
    # Populate band information
    band = 1
    while (True): # Loop until we run out of bands
        bandString = 'Band ' + str(band) + ' Block='
        bandLoc = textOutput.find(bandString)
        if bandLoc < 0: # Ran out of bands
            break
        
        # Found the band, read pertinent information
        bandInfo = {}
        
        # Get the type string
        bandLine = IrgStringFunctions.getLineAfterText(textOutput, bandString)
        typePos  = bandLine.find('Type=')
        commaPos = bandLine.find(',')
        typeName = bandLine[typePos+5:commaPos-1]
        bandInfo['type'] = typeName
        
        outputDict['band_info'] = bandInfo
        
        band = band + 1 # Move on to the next band
        
    
    
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
    cmd = ['geoRefTool --printBounds', geoTiffPath]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    textOutput, err = p.communicate()

    # Check that the call did not fail
    if (textOutput.find('Failed') >= 0):
        raise Exception('Error: getGeoTiffBoundingBox failed on input image: ' + geoTiffPath)
    
    # Parse the output
    minLat = float( IrgStringFunctions.getLineAfterText(textOutput, 'Min latitude  =') )
    maxLat = float( IrgStringFunctions.getLineAfterText(textOutput, 'Max latitude  =') )
    minLon = float( IrgStringFunctions.getLineAfterText(textOutput, 'Min longitude =') )
    maxLon = float( IrgStringFunctions.getLineAfterText(textOutput, 'Max longitude =') )
    
    return (minLon, maxLon, minLat, maxLat)


def getImageBoundingBox(filePath):
    """Returns (minLon, maxLon, minLat, maxLat) for a georeferenced image file"""

    extension = os.path.splitext(filePath)[1]
    if '.cub' in extension:
        return IrgIsisFunctions.getIsisBoundingBox(filePath)
    else: # Handle all other types
        return getGeoTiffBoundingBox(filePath)
          
    # Any other file types will end up raising some sort of exception
    
    
    
    

def build_vrt( fullImageSize, tileLocs, tilePaths, outputPath ):
    """Generates a VRT file from a set of image tiles and their locations in the output image"""

    outputFolder = os.path.dirname(outputPath)

    f = open(outputPath, 'w')
    f.write("<VRTDataset rasterXSize=\"%i\" rasterYSize=\"%i\">\n" % (int(fullImageSize[0]),int(fullImageSize[1])) ) # Write whole image size

    #
    ## If a tile is missing, for example, in the case we
    ## skipped it when it does not intersect user's crop box,
    ## substitute it with a different one, to ensure the mosaic
    ## does not have holes. --> Does this make sense?
    #goodFilename = ""
    #for tile in tiles: # Find the first valid tile (we don't care which one)
    #    directory = settings['out_prefix'][0] + tile.name_str()
    #    filename  = directory + "/" + tile.name_str() + tile_postfix
    #    if os.path.isfile(filename):
    #        goodFilename = filename
    #        break
    #if goodFilename == "":
    #    raise Exception('No tiles were generated')

    
    # Read some metadata from one of the tiles
    gdalInfo = getImageGeoInfo(tilePaths[0])
    
    num_bands = len(gdalInfo['band_info'])
    data_type = gdalInfo['band_info'][0]['type']

    # This special metadata value is only used for ASP stereo point cloud files!    
    if 'point_offset' in gdalInfo:
        f.write("  <Metadata>\n    <MDI key=\"" + 'POINT_OFFSET' + "\">" +
                gdalInfo['point_offset'][0] + "</MDI>\n  </Metadata>\n")
      

    # Write each band
    for b in range( 1, num_bands + 1 ):
        f.write("  <VRTRasterBand dataType=\"%s\" band=\"%i\">\n" % (data_type, b) ) # Write band data type and index

        for tile, tileLoc in zip(tilePaths, tileLocs):
            filename  = tile
            
            imageSize = getImageSize(filename) # Get the image size for this tile

            ## Replace missing tile paths with the good tile we found earlier
            #if not os.path.isfile(filename): filename = goodFilename

            relative = os.path.relpath(filename, outputPath) # Relative path from the output file to the input tile
            f.write("    <SimpleSource>\n")
            f.write("       <SourceFilename relativeToVRT=\"1\">%s</SourceFilename>\n" % relative) # Write relative path
            f.write("       <SourceBand>%i</SourceBand>\n" % b)
            f.write("       <SrcRect xOff=\"%i\" yOff=\"%i\" xSize=\"%i\" ySize=\"%i\"/>\n" % (tileLoc[0], tileLoc[1], imageSize[0], imageSize[1]) ) # Source ROI (entire tile)
            f.write("       <DstRect xOff=\"%i\" yOff=\"%i\" xSize=\"%i\" ySize=\"%i\"/>\n" % (tileLoc[0], tileLoc[1], imageSize[0], imageSize[1]) ) # Output ROI (entire tile)
            f.write("    </SimpleSource>\n")
        f.write("  </VRTRasterBand>\n")
    f.write("</VRTDataset>\n")
    f.close()    
    
    
    
    
    
    
    
    
    
    
    

