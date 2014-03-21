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

import IrgGeoFunctions, IrgIsisFunctions


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Generates a kml overlay for a cube
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def generateKml(boundingBox, outputPath):
    
    # Initialize kml document
    kml = simplekml.Kml()
    kml.document.name = outputPath
    kml.hint = 'target=moon'

    poly = kml.newpolygon(name=outputPath, outerboundaryis=[(boundingBox[0], boundingBox[3]),  # Top left
                                                            (boundingBox[1], boundingBox[3]),  # Top right
                                                            (boundingBox[1], boundingBox[2]),  # Bot right
                                                            (boundingBox[0], boundingBox[2]),  # Bot left
                                                            (boundingBox[0], boundingBox[3])]) # Top right
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
        bb = IrgIsisFunctions.getIsisBoundingBox(cubePath)
        
        # Generate a kml plot of the cube
        generateKml(bb, outputPath)

        print "Finished"
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2


if __name__ == "__main__":
    sys.exit(main())
