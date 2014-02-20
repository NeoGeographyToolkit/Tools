#!/usr/bin/env python
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

import os, glob, optparse, re, shutil, subprocess, sys, string, time, simplekml


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Generates a kml overlay for a cube
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


    

# TODO: Merge this with IsisTools functions
def getCubeSize(cubePath):
    """Returns the size [samples, lines] in a cube"""

    # Make sure the input file exists
    if not os.path.exists(cubePath):
        raise Exception('Cube file ' + cubePath + ' not found!')
       
    # Use subprocess to suppress the command output
    cmd = ['gdalinfo', cubePath]
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

    
# TODO: Merge this with IsisTools functions
def getPixelLocInCube(cubePath, sample, line, workDir=''):
    """Returns the BodyFixedCoordinate of a pixel from a cube"""

    # Make sure the input file exists
    if not os.path.exists(cubePath):
        raise Exception('Cube file ' + cubePath + ' not found!')

    # Default working directory is the cubePath folder
    outputFolder = workDir
    if workDir == '':
        outputFolder = os.path.dirname(cubePath)
       
    if not (len(outputFolder) > 1) and os.path.exists(outputFolder):
        os.mkdir(outputFolder)

    # Call ISIS campt function to compute the pixel location
    tempTextPath = os.path.join(outputFolder, 'camptOutput.txt')
    if os.path.exists(tempTextPath):
        os.remove(tempTextPath) # Make sure any existing file is removed!
        
    # Use subprocess to suppress the command output
    cmd = ['campt', 'from=', cubePath, 'to=', tempTextPath, 'sample=', str(sample), 'line=', str(line)]
    FNULL = open(os.devnull, 'w')
    subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT)

    # Check that we created the temporary file
    if not os.path.exists(tempTextPath):
        raise Exception('campt failed to create temporary file ' + tempTextPath)
    
    infoFile = open(tempTextPath, 'r')
        
    # Read in the output file to extract the pixel coordinates
    gccLine       = ''
    latLine       = ''
    lonLine       = ''
    radiusLine    = ''
    lineAfterBody = False
    for line in infoFile:
        
        # GCC stuff
        if lineAfterBody: # BodyFixedCoordinate takes up two lines
            gccLine       = gccLine + line
            lineAfterBody = False
            
        if (gccLine == ''): # Look for start of the info (this check must come second)
            if (line.find('BodyFixedCoordinate') >= 0):
                gccLine     = line
                lineAfterBody = True
        
        # GDC stuff
        if line.find('PlanetocentricLatitude') >= 0:
            latLine = line
            #print line
        if line.find('PositiveEast180Longitude') >= 0:
            lonLine = line
            #print line
        if line.find('LocalRadius') >= 0:
            radiusLine = line
            #print line

    os.remove(tempTextPath) # Remove the file to clean up

    # Make sure we found the desired lines
    if (gccLine == ''):
        raise Exception("Unable to find BodyFixedCoordinate in file " + tempTextPath)
    if (latLine == ''):
        raise Exception("Unable to find PlanetocentricLatitude in file " + tempTextPath)
    if (lonLine == ''):
        raise Exception("Unable to find PositiveEast180Longitude in file " + tempTextPath)
    if (radiusLine == ''):
        raise Exception("Unable to find LocalRadius in file " + tempTextPath)

    # Extract GCC coordinates
    startParen = gccLine.find('(')
    stopParen  = gccLine.find(')')
    numString  = gccLine[startParen+1:stopParen]
    x,y,z = numString.split(',')

    # Convert output from kilometers to meters
    pixelLocationGcc = [0, 0, 0]
    pixelLocationGcc[0] = float(x) * 1000.0
    pixelLocationGcc[1] = float(y) * 1000.0
    pixelLocationGcc[2] = float(z) * 1000.0
    
    # Extract GDC coordinates
    latStart     = latLine.find('=')+2
    lonStart     = lonLine.find('=')+2
    radiusStart  = radiusLine.find('=')+2
    radiusEnd    = radiusLine.find('<') - 1
    latNumStr    = latLine[latStart:]
    lonNumStr    = lonLine[lonStart:]
    radiusNumStr = radiusLine[radiusStart:radiusEnd]
    pixelLocationGdc = [float(lonNumStr), float(latNumStr), float(radiusNumStr)]
                        
    pixelInformation = dict()
    pixelInformation['gcc'] = pixelLocationGcc
    pixelInformation['gdc'] = pixelLocationGdc
    return pixelInformation

#TODO: Move this into a general location!
def getCubeBoundingBox(cubePath, workDir):
    
    # Get the cube size, then request the positions of the four corners
    cubeSize = getCubeSize(cubePath)
    
    # Note that the underlying ISIS tool is one-based
    topLeft  = getPixelLocInCube(cubePath, 1,           1,           workDir) 
    topRight = getPixelLocInCube(cubePath, cubeSize[0], 1,           workDir)
    botLeft  = getPixelLocInCube(cubePath, 1,           cubeSize[1], workDir)
    botRight = getPixelLocInCube(cubePath, cubeSize[0], cubeSize[1], workDir)
    
    boundingBox = {'topLeft': topLeft['gdc'], 'topRight': topRight['gdc'],
                   'botLeft': botLeft['gdc'], 'botRight': botRight['gdc']}
    
    return boundingBox

def generateKml(boundingBox, outputPath):
    
    # Initialize kml document
    kml = simplekml.Kml()
    kml.document.name = 'test'
    kml.hint = 'target=moon'
        
    poly = kml.newpolygon(name=outputPath, outerboundaryis=[(boundingBox['topLeft' ][0], boundingBox['topLeft' ][1]),
                                                            (boundingBox['topRight'][0], boundingBox['topRight'][1]),
                                                            (boundingBox['botRight'][0], boundingBox['botRight'][1]),
                                                            (boundingBox['botLeft' ][0], boundingBox['botLeft' ][1]),
                                                            (boundingBox['topLeft' ][0], boundingBox['topLeft' ][1])])
#    poly.style.polystyle.color   = simplekml.Color.red # Why is polygon not working?
#        poly.style.polystyle.fill    = 1
#        poly.style.polystyle.outline = 1
    poly.style.linestyle.width   = 5
    
    
    # Save kml document
    kml.save(outputPath)
    
    

def main():

    outputPath = ''

    try:
        try:
            usage = "usage: cube2kml.py [--help][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("-o", "--output-path", dest="outputPath",
                              help="Output path (default replace extension with .kml")
            (options, args) = parser.parse_args()

            if not args: parser.error("need input cube file")
            
            cubePath = args[0]

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        # Determine the output path
        if outputPath == '':
            outputPath = os.path.splitext(cubePath)[0] + '.kml' # Default output path
        outputFolder = os.path.dirname(outputPath)

        # Get the four corners of the cube
        bb = getCubeBoundingBox(cubePath, outputFolder)
        
        # Generate a kml plot of the cube
        generateKml(bb, outputPath)

        print "Finished"
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2


if __name__ == "__main__":
    sys.exit(main())
