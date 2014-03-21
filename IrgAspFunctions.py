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


"""IrgAspFunctions.py - Functions for working with Ames Stereo Pipeline"""

import sys, os, re, subprocess, string, time, errno



def removeIntermediateStereoFiles(stereoPrefix):
    """Deletes intermediate files from a stereo output directory"""

    # List of all non-final-output files
    fileList = ['-align-L.exr', \
                '-align-R.exr', \
                '-DEM.tif.aux.xml', \
                '-D_sub.tif', \
                '-D.tif', \
                '-F.tif', \
                '-GoodPixelMap.tif', \
                '-lMask_sub.tif', \
                '-lMask.tif', \
                '-L_sub.tif', \
                '-L.tif', \
                '-RD.tif' ,\
                '-rMask_sub.tif', \
                '-rMask.tif', \
                '-R_sub.tif' ,\
                '-R.tif']
    # Remove each of those files
    for f in fileList:
        path = stereoPrefix + f
        removeIfExists(path)



def readJitregFile(filePath):
    """Reads the output file from a lronacjitreg call and returns [meanSampleOffset, meanLineOffset]"""

    # Fail if the input file is not present
    if not os.path.isfile(filePath):
        raise Exception('File ' + filePath + ' is missing!')

    averages = [0.0, 0.0]

    f = open(filePath,'r')
    for line in f:
        if ( line.rfind("Average Sample Offset:") >= 0 ):
            index       = line.rfind("Offset:");
            index_e     = line.rfind("StdDev:");
            crop        = line[index+7:index_e];
            if crop == " NULL ": # Check for null value
                raise Exception('Null sample offset in file ' + flat)
            averages[0] = float(crop);
        elif ( line.rfind("Average Line Offset:") >= 0 ):
            index       = line.rfind("Offset:");
            index_e     = line.rfind("StdDev:");
            crop        = line[index+7:index_e];
            if crop == "   NULL ": # Check for null value
                raise Exception('Null sample offset in file ' + flat)
            averages[1] = float(crop);
        elif ( line.rfind("Using IpFind result only:") >= 0 ):
            index       = line.rfind("only:");
            if (line[index + 7] == 1):
                print "Warning: This result based only on IpFind search."
    print str(averages)
    return averages



def getStereoGoodPixelPercentage(inputPrefix, workDir=''):
    """Returns the percentage of good pixels in a stereo output"""

    # Set up input folder
    inputFolder = os.path.dirname(inputPrefix)
    if not os.path.exists(inputFolder):
        raise Exception('Input folder ' + inputFolder + ' not found!')    
    if workDir == '':
        workDir = inputFolder

    
    #TODO: Look for goodPixelMap file!
    
    #TODO: Look for later stage estimates!
    
    # If the later stage files were not found, use the integer correlation file 
    
    # Extract the third band of the D_sub.tif image which contains a good pixel map
    inputPath = inputPrefix + '-D_sub.tif'
    if not os.path.exists(inputPath):
        raise Exception('Could not find file ' + inputPath)
    convertedImagePath = os.path.join(workDir,     'goodPixelMap-D_sub.tif')
    cmd = 'gdal_translate -of GTiff -ot BYTE -b 3 ' + inputPath + ' ' + convertedImagePath
    print cmd
    os.system(cmd)
    
    # Determine the percentage of good pixels   
    cmd = ['gdalinfo', '-hist', convertedImagePath]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    translateOut, err = p.communicate()

    # Parse the gdalinfo output
    bucket  = translateOut.find('buckets')
    colon   = translateOut.find(':', bucket)
    start   = translateOut.find('\n', colon)
    end     = translateOut.find('\n', start+1)
    buckets = translateOut[start+1:end] # Pick off the string containing the buckets
    numbers = buckets.strip().split(' ')
    
    numBad      = int(numbers[0]) # All pixels are in the first (bad) or last (good) buckets
    numGood     = int(numbers[-1])
    percentGood = float(numGood) / float(numGood + numBad)

    return percentGood
  

