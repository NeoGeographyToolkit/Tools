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
    aRadEnd = outputText.find(' ', aRadPos)
    bRadEnd = outputText.find(' ', bRadPos)
    aRad    = float(outputText[aRadPos:aRadEnd])
    bRad    = float(outputText[bRadPos:bRadEnd])
    if aRad != bRad:
        raise Exception('Only spherical ellipsoids are supported!')
    
    # Extract the projection space corner coordinates
    ulPos = outputText.find('Upper Left')
    lrPos = outputText.find('Lower Right')

    ulCoordStart = outputText.find('(', ulPos)
    ulCoordStop  = outputText.find(')', ulPos)
    lrCoordStart = outputText.find('(', lrPos)
    lrCoordStop  = outputText.find(')', lrPos)

    ulCoord = outputText[ulCoordStart+1:ulCoordStop-1]
    lrCoord = outputText[lrCoordStart+1:lrCoordStop-1]

    ulx = float(ulCoord.split(',')[0]) # These are X and Y projection space coordinates
    uly = float(ulCoord.split(',')[1])
    lrx = float(lrCoord.split(',')[0])
    lry = float(lrCoord.split(',')[1])
    
    # Compute distances in projection space
    PI     = 3.14159265358979323846264338327950288419716939937510
    radius = aRad 
    onePiR = radius * PI #5458203.07635 # 180 degree distance along planet's equator
    twoPiR = radius * 2.0 * PI #10916406.153 # Circumference of the Planet's equator


    # Extract the current proj4 string
    proj4Pos       = outputText.find("+proj=")
    proj4End       = outputText.find("'", proj4Pos) - 1 # Relies on proj4 string ending with a ' symbol
    oldProj4String = outputText[proj4Pos:proj4End]

    # Get the value of +lon_0
    lon0Pos       = oldProj4String.find('+lon_0=') + 7
    lon0End       = oldProj4String.find(' ', lon0Pos)
    oldLon0String = oldProj4String[lon0Pos-7:lon0End]
    lon0          = float(oldProj4String[lon0Pos:lon0End])

    if lon0 == 0.0: # No need to change the proj4 string
        newProj4String = oldProj4String
    else:
        # Need to change the lon_0 value to zero and account for it
        newProj4String = oldProj4String.replace(oldLon0String, '+lon_0=0')
        
        # This amount is being implicitly added to the projected coordinates
        lon0Correction = lon0 * (PI / 180) * radius
        ulx += lon0Correction
        lrx += lon0Correction


    # In order for many programs to ingest data, the projected space coordinates must
    #  be in the +/- 180 degree (onePiR) range.

    if (abs(ulx) < onePiR) and (abs(lrx) < onePiR):
        print 'Both images are in the +/-180 degree range, no correction needed!'
        return True

    if (abs(ulx) < onePiR) or (abs(lrx) < onePiR):
        raise Exception('How to handle an image that straddles 180???')

    # The entire image is outside +/-180, now test which way
    if ulx < onePiR:
        offset = twoPiR
    else:
        offset = twoPiR * -1

    # Perform an in-place correction on the file
    cmd = '/home/pirl/smcmich1/programs/gdal-1.11.0-install/bin/gdal_edit.sh ' + tiffPath + ' -a_srs "'+newProj4String +'" -a_ullr '+ str(ulx+offset) +' '+ str(uly) +' '+ str(lrx + offset) +' '+ str(lry)
    #cmd = 'gdal_edit.py ' + tiffPath + ' -a_srs "'+newProj4String +'" -a_ullr '+ str(ulx+offset) +' '+ str(uly) +' '+ str(lrx + offset) +' '+ str(lry)
    print cmd
    os.system(cmd)

    # We could use gdal_translate to do this but it would not be in-place.  This is an issue
    #   for large images we don't want to copy.

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

        if (options.normEqcLon):
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
