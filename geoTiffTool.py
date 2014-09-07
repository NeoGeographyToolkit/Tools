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

import IrgGeoFunctions, IrgIsisFunctions, IrgSystemFunctions


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Interface for working with GeoTiff files
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg



def correctTifPast180(tiffPath):
    '''Normalizes the projection space location of an image beyond 180 degrees'''

    # Todo: Check for gdal_edit.py

    # First figure out if the file is past 180 degrees
    cmd = ['gdalinfo',  tiffPath, '-proj4']
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    outputText, err = p.communicate()

    # Verify projection type
    if not '+proj=eqc' in outputText:
        raise Exception('Only EQC projection files are supported!')

    # Get the planetary radius out of the proj4 string
    aRadPos = outputText.find('+a=') + 3
    bRadPos = outputText.find('+b=') + 3
    aRadEnd = outputText.find(' ', arRadPos)
    bRadEnd = outputText.find(' ', brRadPos)
    aRad    = float(outputText[aRadPos:aRadEnd])
    bRad    = float(outputText[bRadPos:bRadEnd])
    if aRad != bRad:
        raise Exception('Only spherical ellipsoids are supported!')
    

    ulPos = outputText.find('Upper Left')
    lrPos = outputText.find('Lower Right')

    ulCoordStart = outputText.find('(', ulPos)
    ulCoordStop  = outputText.find(')', ulPos)
    lrCoordStart = outputText.find('(', lrPos)
    lrCoordStop  = outputText.find(')', lrPos)

    ulCoord = outputText[ulCoordStart+1:ulCoordStop-1]
    lrCoord = outputText[lrCoordStart+1:lrCoordStop-1]

    ulx = float(ulCoord.split(',')[0])
    uly = float(ulCoord.split(',')[1])
    lrx = float(lrCoord.split(',')[0])
    lry = float(lrCoord.split(',')[1])
    
    PI     = 3.14159265358979323846264338327950288419716939937510
    radius = aRad 
    onePiR = radius * PI #5458203.07635 # 180 degree distance along planet's equator
    twoPiR = radius * two * PI #10916406.153 # Circumference of the Planet's equator

    if (abs(ulx) < onePiR) and (abs(lrx) < onePiR):
        # Both images are in the +/-180 degree range, no correction needed!
        return True

    if (abs(ulx) < onePiR) or (abs(lrx) < onePiR):
        raise Exception('How to handle an image that straddles 180???')

    # The entire image is outside +/-180, now test which way
    if ulx < onePiR:
        offset = twoPiR
    else:
        offset = twoPiR * -1

    # Perform an in-place correction on the file
    #cmd = '~/programs/gdal_install/bin/gdal_edit.sh ' + tiffPath +' -a_ullr '+ str(ulx+offset) +' '+ str(uly) +' '+ str(lrx + offset) +' '+ str(lry)
    cmd = 'gdal_edit.py ' + tiffPath +' -a_ullr '+ str(ulx+offset) +' '+ str(uly) +' '+ str(lrx + offset) +' '+ str(lry)
    print cmd
    os.system(cmd)

    return 0
    

def main():

    outputPath = ''

    try:
        try:
            usage = "usage: geoTiffTool.py [--help][--manual] geotiffPath\n"# [lat] [lon] [lat2] [lon2]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
  
            ## Many GeoTiff images are stored in projected coordinates which makes this tricky
            #parser.add_option("--crop", action="store_true", dest="doCrop", default=False,
            #                            help="Crops the image between the four lat/lon arguments.")

            parser.add_option("--normalize-eqc-lon", action="store_true", dest="normEqcLon", default=False,
                                        help="Normalizes an EQC image to lie between +/- 180 degrees.")

            # Add other features as the need arises

            #parser.add_option("-o", "--output-path", dest="outputPath",
            #                        help="Output path (default replace extension with .kml")
            
            (options, args) = parser.parse_args()

            # Parse positional arguments
            if not args: parser.error("need input cube file")           
            tiffPath = args[0]
            #if len(args) > 1:
            #    lat = float(args[1])
            #if len(args) > 2:
            #    lon = float(args[2])
            #if len(args) > 3:
            #    lat2 = float(args[3])
            #if len(args) > 4:
            #    lon2 = float(args[4])

        except optparse.OptionError, msg:
            raise Usage(msg)

        if (normEqcLon):
            return correctTifPast180(tiffPath)


        #if doCrop:
        #    if len(args) > 5:
        #        raise Exception('Not enough arguments supplied!  Need input + two lat/lon pairs')
        #    
        #    # TODO: Verify that the geotiff image is not already in GDC coordinates!
        #    
        #    # Get the pixel coordinates corresponding to the lat/lon coordinates
        #    cmd = ['geoRefTool', '--lat', lat, '--lon', lon, '--outputFormat' '2']
        #    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        #    projectionInfo, err = p.communicate()

  
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2


if __name__ == "__main__":
    sys.exit(main())
