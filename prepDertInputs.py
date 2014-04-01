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

import sys, os, glob, optparse, re, shutil, subprocess, string, time, logging, threading, math


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Generate a DERT landscape directory from orthoimage and DEM inputs
'''
    sys.exit()


# TODO: Move these functions to a common location!

def readPairInParen(string, startPos):
    """Reads a pair of numbers contained in parenthesis like this: (3413.55, 4103.456)"""
    
    # Find the bounds
    startParen = string.find('(', startPos)
    commaLoc   = string.find(',', startParen)
    stopParen  = string.find(')', commaLoc)
    
    # Extract the numbers
    num1 = float(string[startParen+1:commaLoc-1])
    num2 = float(string[commaLoc+1  :stopParen-1])
    
    return (num1, num2)

def getGdalBounds(filePath):

    # Make sure the input file exists
    if not os.path.exists(filePath):
        raise Exception('Error: file ' + filePath + ' not found!')

    # Use subprocess to read the command output
    cmd = ['gdalinfo', filePath]
    FNULL = open(os.devnull, 'w')
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    cmdOut, err = p.communicate()
    
    ulStart = cmdOut.find('Upper Left' )
    llStart = cmdOut.find('Lower Left' )
    urStart = cmdOut.find('Upper Right')
    lrStart = cmdOut.find('Lower Right')
    
    ulVals = readPairInParen(cmdOut, ulStart)
    llVals = readPairInParen(cmdOut, llStart)
    urVals = readPairInParen(cmdOut, urStart)
    lrVals = readPairInParen(cmdOut, lrStart)
    
    # Return minX, minY, maxX, maxY
    return (ulVals[0], lrVals[1], lrVals[0], ulVals[1])

#TODO: Really need a python bounding box class
def getBoundsOverlap(bb1, bb2):
    """Returns the intersection of two bounding boxes"""

    minX = max(bb1[0], bb2[0])
    minY = max(bb1[1], bb2[1])
    maxX = min(bb1[2], bb2[2])
    maxY = min(bb1[3], bb2[3])

    return (minX, minY, maxX, maxY)

def cropToGeoBounds(inputPath, outputPath, bounds):
    """Use gdal_translate to crop a file to the given bounds"""
    
    if os.path.exists(outputPath):
        return True
    
    cmd = ('gdal_translate ' + inputPath + ' ' + outputPath +
                            ' -projwin ' + str(bounds[0]) + ' ' + str(bounds[3]) +
                                     ' ' + str(bounds[2]) + ' ' + str(bounds[1]) )
    print cmd
    os.system(cmd)
    
    return os.path.exists(outputPath)

def getImagePixelSize(inputPath):
    """Use gdalinfo to get the pixel size"""
    
    # Use subprocess to read the command output
    cmd = ['gdalinfo', inputPath]
    FNULL = open(os.devnull, 'w')
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    cmdOut, err = p.communicate()
    
    lineStart = cmdOut.find('Pixel Size')
    eqPos     = cmdOut.find('=', lineStart)
        
    pixelSize = readPairInParen(cmdOut, eqPos)
    return pixelSize[0] # Currently assume x and y size are equal!

def getLowerPowerOfTwo(number):
    """Returns the next power of two lower than the given number"""
    return int(math.log(number, 2))

def computeDemRescale(demPixelSize, imagePixelSize):
    """Compute required DERT rescale of the DEM\n
       Size of the image must be (power of 2) times the DEM size"""
       
    startRatio      = demPixelSize / imagePixelSize
    
    lowPowerOfTwo   = getLowerPowerOfTwo(startRatio)*2
    
    #TODO: Handle case where the image is lower res than the DEM!
   
    print 'Output image is ' + str(lowPowerOfTwo) + ' times the size of the output DEM'
    
    # Note: larger pixel size here means lower number of image pixels
    newDemPixelSize = imagePixelSize * lowPowerOfTwo
    
    newDemPixelSizeFraction = newDemPixelSize / demPixelSize
    
    return newDemPixelSizeFraction
    
# TODO: This is a duplicate of code in IsisTools.py!
def getImageSize(imagePath):
    """Returns the size [samples, lines] in a cube"""

    # Make sure the input file exists
    if not os.path.exists(imagePath):
        raise Exception('Image file ' + imagePath + ' not found!')
       
    # Use subprocess to suppress the command output
    cmd = ['gdalinfo', imagePath]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    textOutput, err = p.communicate()

    # Extract the size from the text
    sizePos    = textOutput.find('Size is')
    endPos     = textOutput.find('\n', sizePos+7)
    sizeStr    = textOutput[sizePos+7:endPos]
    sizeStrs   = sizeStr.strip().split(',')
    numSamples = int(sizeStrs[0])
    numLines   = int(sizeStrs[1])
    
    size = [numSamples, numLines]
    return size

    
def resizeImage(inputPath, outputPath, ratio):
    """Resizes an image on disk by the given ratio"""
    
    if os.path.exists(outputPath):
        return True
    
    inputImageSize = getImageSize(inputPath)
    
    newSize = (inputImageSize[0] * ratio, inputImageSize[1] * ratio)
    
    cmd = ('gdalwarp -r bilinear -ts ' + str(newSize[0]) + ' ' + str(newSize[1]) +
                     ' ' + inputPath + ' ' + outputPath)
    print cmd
    os.system(cmd)
    
    return os.path.exists(outputPath)

def callLandscapefactory(imagePath, demPath, outputFolder):

    # TODO: Do we need to manually handle nodata values?

    # Prep the DEM portion    
    cmd = ('landscapefactory -landscape=' + outputFolder +
                           ' -file=' + demPath + ' -layer=elevation -tilesize=128 -globe=moon')
    print cmd
    os.system(cmd)
    
    print '\n------------------------\n'
    
    # Prep the image portion    
    cmd = ('landscapefactory -landscape=' + outputFolder +
                           ' -file=' + imagePath + ' -layer=grayphoto -tilesize=128 -globe=moon')
    print cmd
    os.system(cmd)
    
    return True

# Set up a pair of files to be loaded into DERT
def main(argsIn):
    
    try:
        usage = "usage: prepDertInputs.py [--help][--manual]\n  "
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("--manual", action="callback", callback=man,
                          help="Read the manual.")
        
        parser.add_option("--dem",          dest="demPath",      help="Path to input DEM")
        parser.add_option("--image",        dest="imagePath",    help="Path to input orthoimage")
        parser.add_option("--outputFolder", dest="outputFolder", help="Output folder to write to")
        
#        parser.add_option("--force", action="store_true", dest="force",
#                          help="Force overwrite of existing files.")
        
        parser.add_option("--keep", action="store_true", dest="keep",
                          help="Retain intermediate files.")
            
        (options, args) = parser.parse_args()

        if not options.demPath or not options.imagePath or not options.outputFolder:
            parser.error('Missing required inputs!')

    except optparse.OptionError, msg:
        raise Usage(msg)

    # Generate temporary file paths
    convertedDemPath   = os.path.join(options.outputFolder, 'tempDertDem.tif')
    convertedImagePath = os.path.join(options.outputFolder, 'imageForDert.tif')   
    resizedDemPath     = os.path.join(options.outputFolder, 'demForDert.tif')
        
    if not os.path.exists(options.outputFolder):
        os.mkdir(options.outputFolder)
    
    print 'Starting prepDertInputs.py'
    
    #TODO: Verify they are in the same coordinate system!
    
    print 'Determining common bounding box...'
    
    # The projected image and the DEM need to cover the exact same area
    bbDem    = getGdalBounds(options.demPath)
    bbImage  = getGdalBounds(options.imagePath)
    bbCommon = getBoundsOverlap(bbDem, bbImage)
    
    print 'Found common bounding box: ' + str(bbCommon)
    
    print 'Cropping inputs to common bounding box...'
    
    cropToGeoBounds(options.demPath,   convertedDemPath,   bbCommon)
    cropToGeoBounds(options.imagePath, convertedImagePath, bbCommon)
    
    print 'Computing required DEM rescale...'
    
    # The image needs to be [power of two] times higher resolution than the DEM
    # - We shrink the DEM down to achieve this
    demPixelSize       = getImagePixelSize(convertedDemPath)
    imagePixelSize     = getImagePixelSize(convertedImagePath)
    requiredDemRescale = computeDemRescale(demPixelSize, imagePixelSize)
    
    print 'Rescaling DEM...'
    
    resizeImage(convertedDemPath, resizedDemPath, requiredDemRescale)

    print 'Calling landscapefactory...'
    callLandscapefactory(convertedImagePath, resizedDemPath, options.outputFolder)

    # Clean up
    if not options.keep:
        os.remove(convertedDemPath)
        os.remove(convertedImagePath)
        os.remove(resizedDemPath)

    print 'Finished prepDertInputs.py!'


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
    