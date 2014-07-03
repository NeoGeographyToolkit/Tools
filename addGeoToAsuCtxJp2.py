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

import sys, os, glob, optparse, re, shutil, subprocess, string, time


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Merge the geo information with an ASU map projected CTX JP2 file.
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#--------------------------------------------------------------------------------

# TODO: Move to IrgFileFunctions!
def replaceFileExtension(inputPath, ext):
    '''Replace the extension on a file'''
    return os.path.splitext(inputPath)[0] + ext
    

def addGeoData(inputJp2Path, inputHeaderPath, outputPath, keep=False):
    """Does the actual work of adding the geo data"""

    if not os.path.exists(inputJp2Path):
        raise Exception('Input file ' + inputJp2Path + ' does not exist!')
    if not os.path.exists(inputHeaderPath):
        raise Exception('Input file ' + inputHeaderPath + ' does not exist!')


    # Get the needed paths
    prjPath = replaceFileExtension(inputJp2Path, '.prj')
    vrtPath = replaceFileExtension(inputJp2Path, '.vrt')
    
    # The perl script works best from the input folder
    originalDirectory = os.getcwd()
    print originalDirectory
    inputFolder = os.path.dirname(inputJp2Path)
    print inputFolder
    if inputFolder != '':
        os.chdir(inputFolder)

    # Call perl script, then return to original directory
    cmd = 'isis3world.pl -J -prj ' + os.path.basename(inputHeaderPath)
    print cmd
    os.system(cmd)
    if inputFolder != '':
        os.chdir(originalDirectory)

    # If the first command is not called from the input folder, need to do this
    #mv .j2w <INPUT FOLDER>/B01_009838_2108_XI_30N319W.j2w
    #mv .prj <INPUT FOLDER>/B01_009838_2108_XI_30N319W.prj

    if (not os.path.exists(prjPath)):
        raise Exception('isis3world.pl script failed to create required output files!')

    # Next conversion command
    cmd = 'gdal_translate -of VRT -a_srs ESRI::"'+ prjPath +'" -a_nodata 0  '+ inputJp2Path +' '+ vrtPath
    print(cmd)
    os.system(cmd)

    # Finish the conversion
    cmd = 'gdal_translate -of JP2OpenJPEG '+ vrtPath +' '+ outputPath
    print(cmd)
    os.system(cmd)

    # Clean up temporary files
    if not keep:
        os.remove(vrtPath)
        os.remove(prjPath)

    return True

def main():

    print '#################################################################################'
    print "Running addGeoToAsuCtxJp2.py"

    try:
        try:
            usage = "usage: addGeoToAsuCtxJp.py <inputPath> <outputPath> [--keep][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)

            parser.set_defaults(keep=False)

            parser.add_option("--label", dest="inputHeaderFile", default="",
                              help="Path to the label file.")

            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("--keep", action="store_true", dest="keep",
                              help="Do not delete the temporary files.")
            (options, args) = parser.parse_args()
            
            if len(args) < 1:
                parser.error('Missing required input!')
            options.inputJp2Path = args[0]
            if len(args) > 1:
                options.outputPath = args[1]
            else:
                options.outputPath = replaceFileExtension(options.inputJp2Path, '.tif')
            
            # If the path to the header file was not provided, assume the default naming convention
            if options.inputHeaderFile == "":
                options.inputHeaderFile = replaceFileExtension(options.inputJp2Path, '.scyl.isis.hdr')
            
        except optparse.OptionError, msg:
            raise Usage(msg)

        startTime = time.time()

        addGeoData(options.inputJp2Path, options.inputHeaderFile, options.outputPath, options.keep)

        endTime = time.time()

        print "Finished in " + str(endTime - startTime) + " seconds."
        print '#################################################################################'
        return 0

    except Usage, err:
        print err
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main())
