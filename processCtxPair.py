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
Generates a stereo DEM from two LRONAC pairs using SBA and LOLA for increased accuracy.
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#--------------------------------------------------------------------------------

def replaceExtensionAndFolder(inputPath, outputFolder, newExtension):
    nexExt = os.path.splitext(inputPath)[0] + newExtension
    return os.path.join(outputFolder, os.path.basename(newExt))


def prepareImage(inputPath, workDir, keep):
    """Prepare a single CTX image for processing"""

    # Set up paths
    cubPath = replaceExtensionAndFolder(inputPath, workDir, '.cub')
    calPath = replaceExtensionAndFolder(inputPath, workDir, '.cal.cub')

    # Convert to ISIS format
    cmd = 'mroctx2isis from=' + inputImagePath  + ' to=' + cubPath
    os.system(cmd)
    
    # Init Spice data
    cmd = 'spiceinit from=' + cubPath
    os.system(cmd)
    
    # Apply image correction
    cmd = 'ctxcal from='+cubPath+' to='+calPath
    os.system(cmd)

    #you can also optionally run} ctxevenodd \textnormal{on the} cal.cub \textnormal{files, if needed}

    if not keep:
        os.remove(cubPath)
    
    return calPath

def main():

    print '#################################################################################'
    print "Running processCtxPair.py"

    try:
        try:
            usage = "usage: processCtxPair.py <left image> <right image> <output prefix> [--workDir <folder>][--keep][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)

            parser.set_defaults(keep=False)

            parser.add_option("--workDir",  dest="workDir",  help="Folder to place intermediate files in")

            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("--keep", action="store_true", dest="keep",
                              help="Do not delete the temporary files.")
            (options, args) = parser.parse_args()
            
            if len(args) < 3:
                parser.error('Missing required input!')
            options.leftImagePath  = args[0]
            options.rightImagePath = args[1]
            options.outputPrefix   = args[2]
            
            if not options.workDir:
                options.workDir = os.path.dirname(options.outputPath)

        except optparse.OptionError, msg:
            raise Usage(msg)

        startTime = time.time()

        # Do individual input image preparations
        leftCalPath  = prepareImage(options.leftPath,  options.workDir, options.keep)
        rightCalPath = prepareImage(options.rightPath, options.workDir, options.keep)
        
        # Do joint prepration
        cmd = 'cam2map4stereo.py ' + leftCalPath + ' ' + rightCalPath
        os.system(cmd)
        leftMapPath  = replaceExtensionAndFolder(options.leftPath,  workDir, '.map.cub')
        rightMapPath = replaceExtensionAndFolder(options.rightPath, workDir, '.map.cub')
  
        # Final stereo call
        cmd = ('parallel_stereo.py ' + leftMapPath + ' ' + rightMapPath + ' ' + options.outputPrefix
               + ' --alignment affineepipolar --subpixel-mode 3 --corr-timeout 400'
               + ' --filter-mode 1 --subpixel-max-levels 0')
        os.system(cmd)

        # Clean up temporary files
        if not options.keep:
            os.remove(leftCalPath)
            os.remove(rightCalPath)
            os.remove(leftMapPath)
            os.remove(rightMapPath)

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
