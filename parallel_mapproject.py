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

import sys

import os, glob, re, shutil, subprocess, string, time, errno, optparse, math

import IrgFileFunctions, IrgIsisFunctions, IrgPbsFunctions, IrgSystemFunctions

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Calls mapproject in parallel with ISIS cameras.
'''
    sys.exit()


def generateTileName(startX, startY, stopX, stopY):
    """Generate the name of a tile based on its location"""
    
    tileString = 'tile_' + str(startX) + '_' + str(startY) + '_' + str(stopX) + '_' + str(stopY) + '_.tif'
    return tileString
    

def generateTileList(fullWidth, fullHeight, tileSize):
    """Generate a full list of tiles for this image"""

    numTilesX = int(math.ceil(fullWidth  / float(tileSize)))
    numTilesY = int(math.ceil(fullHeight / float(tileSize)))

    tileList = []
    for r in range(0, numTilesY):
        for c in range(0, numTilesX):
            
            # Starting pixel positions for the tile
            tileStartY = r * tileSize
            tileStartX = c * tileSize
            
            # Determine the size of this tile
            thisWidth  = tileSize
            thisHeight = tileSize
            if (r == numTilesY-1): # If the last row
                thisHeight = fullHeight - tileStartY # Height is last remaining pixels
            if (c == numTilesX-1): # If the last col
                thisWidth  = fullWidth  - tileStartX # Width is last remaining pixels
            
            # Get the end pixels for this tile
            tileStopY  = tileStartY + thisHeight # Stop values are exclusive
            tileStopX  = tileStartX + thisWidth
            
            # Create a name for this tile
            # - Tile format is tile_col_row_width_height_.tif
            tileString = generateTileName(tileStartX, tileStartY, tileStopX, tileStopY)
            
            tileList.append((tileStartX, tileStartY, tileStopX, tileStopY, tileString))
    
    return (numTilesX, numTilesY, tileList)

def handleArguments(args):
    """Split up arguments into required and optional lists which will be passed to subprocess"""

    requiredList = []
    optionsList  = []
       
    # Loop through all entries.
    iterable = iter(range(0, len(args)))
    for i in iterable:
        a = args[i]
        if (i < len(args)-1): # Don't load the next value when we are at the end!
            n = args[i+1]
        else:
            n = '-' # This will just cause us to get out of the loop
        
        if IrgSystemFunctions.isCmdOption(a):      # This is the start of an option.
            optionsList.append(a)  # Record this entry.
            
            if IrgSystemFunctions.isCmdOption(n):  # The next entry is the start of another option so this one has no values.
                continue
            
            optionsList.append(n)  # Otherwise record the next entry as a value.
            iterable.next()              # Skip the next entry in the loop.

            if (a == '--t_projwin') or (a == '--t_pixelwin'):  # These arguments have four values, not just one.
                optionsList.append(args[i+2])              # Add the additional three arguments and skip them in the loop.
                optionsList.append(args[i+3])
                optionsList.append(args[i+4])
                iterable.next()
                iterable.next()
                iterable.next()
        
        else: # This is one of the three positional arguments
            requiredList.append(a)
    
    # Return the two lists
    return (requiredList, optionsList)

def writeSingleTile(options):
    """Writes a single tile according to the options"""

    # Determine the name of the tile we need to write
    tileName = generateTileName(options.pixelStartX, options.pixelStartY, options.pixelStopX, options.pixelStopY)
    tilePath = os.path.join(options.workDir, tileName)
       
    # Just call the command for a single tile!
    cmd = ['mapproject',  '--t_pixelwin', str(options.pixelStartX), str(options.pixelStartY), str(options.pixelStopX), str(options.pixelStopY),
                               options.demPath, options.imagePath, tilePath]
    cmd = cmd + options.extraArgs # Append other options
    IrgSystemFunctions.executeCommand(cmd, options.suppressOutput)
      
    if options.convertTiles: # Make uint8 version of the tile for debugging
        
        tilePathU8 = os.path.splitext(tilePath)[0] + 'U8.tif'
        cmd = ['gdal_translate', '-ot', 'byte', '-scale', tilePath, tilePathU8]
        IrgSystemFunctions.executeCommand(cmd, options.suppressOutput)

    return 0

#------------------------------------------------------------------------------

def main(argsIn):

    try:
        usage      = "usage: parallel_mapproject.py [options] <dem> <camera-image> <output>"
        epilogText = "This also accepts all 'mapproject' arguments though the 'threads' argument will typically be ignored."
        parser     = IrgSystemFunctions.PassThroughOptionParser(usage=usage, epilog=epilogText) # Use parser that ignores unknown options


        parser.add_option("--num-processes",  dest="numProcesses", type='int', default=None,
                                              help="Number of processes to use (default program tries to choose best)")

        parser.add_option('--nodes-list',  dest='nodesListPath', default=None,
                                           help='The list of computing nodes, one per line. ' + \
                                                'If not provided, run on the local machine.')

        parser.add_option('--tile-size',  dest='tileSize', default=1000,
                                           help='Size of square tiles to break up processing in to.')


        # Directory where the job is running
        parser.add_option('--work-dir',  dest='workDir', default=None,
                                         help='Working directory to assemble the tiles in')

        parser.add_option("--suppress-output", action="store_true", default=False,
                                               dest="suppressOutput",  help="Suppress output of sub-calls.")

        parser.add_option("--manual", action="callback", callback=man,
                                       help="Read the manual.")
               
        # DEBUG options
        parser.add_option("--keep", action="store_true", dest="keep", default=False,
                                    help="Do not delete the temporary files.")
        parser.add_option("--convert-tiles",  action="store_true", dest="convertTiles",
                                              help="Generate a uint8 version of each tile")



        ## Debug options
        #p.add_option('--dry-run',   dest='dryrun', default=False, action='store_true',
        #                            help=optparse.SUPPRESS_HELP)
        #p.add_option('--verbose',   dest='verbose', default=False, action='store_true',
        #                            help=optparse.SUPPRESS_HELP)        
    
    
        # PRIVATE options
        # These specify the tile location to request, bypassing the need to query mapproject.
        parser.add_option('--pixelStartX', dest='pixelStartX', default=None, type='int',
                                           help=optparse.SUPPRESS_HELP)
        parser.add_option('--pixelStartY', dest='pixelStartY', default=None, type='int',
                                           help=optparse.SUPPRESS_HELP)
        parser.add_option('--pixelStopX',  dest='pixelStopX', default=None, type='int',
                                           help=optparse.SUPPRESS_HELP)
        parser.add_option('--pixelStopY',  dest='pixelStopY', default=None, type='int',
                                           help=optparse.SUPPRESS_HELP)



        # This call handles all the parallel_mapproject specific options.
        (options, args) = parser.parse_args(argsIn)

        # This will parse all the mapproject options.
        requiredList, optionsList = handleArguments(args)

        # Check the required positional arguments.
        if len(requiredList) < 1:
            parser.error("Need path to DEM")
        if len(requiredList) < 2:
            parser.error("Need path to input image")
        if len(requiredList) < 3:
            parser.error("Need output path")

        options.demPath    = requiredList[0]
        options.imagePath  = requiredList[1]
        options.outputPath = requiredList[2]


        # Any additional arguments need to be forwarded to the mapproject function
        options.extraArgs = optionsList

    except optparse.OptionError, msg:
        raise Usage(msg)

    startTime = time.time()
    
    # Determine if this is a main copy or a spawned copy
    spawnedCopy = ( (options.pixelStartX is not None) and (options.pixelStartY is not None) and
                    (options.pixelStopX  is not None) and (options.pixelStopY  is not None) and options.workDir )
    
    if spawnedCopy: # This copy was spawned to process a single tile
        
        return writeSingleTile(options) # Just call a function to handle this and then we are done!   

    # If the input image is NOT an ISIS image AND we are running on a single machine we can
    #  just use the multi-threading capability of the ordinary mapproject call.
    if (not IrgIsisFunctions.isIsisFile(options.imagePath)) and (not options.nodesListPath):
        cmd = ['mapproject',  options.imagePath, options.demPath, options.outputPath]    
        cmd = cmd + options.extraArgs
        subprocess.call(cmd)
        return 0


    # Otherwise this is the original called process and there are multiple steps to go through
    
    # Call mapproject on the input data using subprocess and record output
    cmd = ['mapproject',  '--query-projection', options.demPath, options.imagePath, options.outputPath]
    cmd = cmd + options.extraArgs # Append other options
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    projectionInfo, err = p.communicate()
    if not options.suppressOutput:
        print projectionInfo
    
        
    # Now find the image size in the output
    startPos    = projectionInfo.find('Output image bounding box:')
    widthStart  = projectionInfo.find('width:', startPos)
    heightStart = projectionInfo.find('height:', widthStart)
    heightEnd   = projectionInfo.find(')', heightStart)
    # Extract values
    fullWidth   = int(projectionInfo[widthStart+7  : heightStart-1])
    fullHeight  = int(projectionInfo[heightStart+8 : heightEnd])
    print 'Output image size is ' + str(fullWidth) + ' by ' + str(fullHeight) + ' pixels.'

    # For now we just break up the image into a user-specified tile size (default 1000x1000)
    numTilesX, numTilesY, tileList = generateTileList(fullWidth, fullHeight, options.tileSize)
    numTiles = numTilesX * numTilesY

    print 'Splitting into ' + str(numTilesX) + ' by ' + str(numTilesY) + ' tiles.'

    # Set up output folder
    outputFolder = os.path.dirname(options.outputPath)
    if outputFolder == '':
        outputFolder = './' # Handle calls in same directory
    outputName   = os.path.basename(options.outputPath)
    IrgFileFunctions.createFolder(outputFolder)
    
    # Make a temporary directory to store the tiles
    if options.workDir:
        tempFolder = options.workDir
    else: # No folder provided, create a default one
        tempFolder = os.path.join(outputFolder, outputName.replace('.', '_') + '_tiles/')
    IrgFileFunctions.createFolder(tempFolder)

    
    # Generate a text file that contains the boundaries for each tile
    argumentFilePath = os.path.join(tempFolder, 'argumentList.txt')
    argumentFile     = file(argumentFilePath, 'w')
    for tile in tileList:
        argumentFile.write(str(tile[0]) + '\t' + str(tile[1]) + '\t' + str(tile[2]) + '\t' + str(tile[3]) + '\n')
    argumentFile.close()
    
    # Indicate to GNU Parallel that there are multiple tab-seperated variables in the text file we just wrote
    parallelArgs = ['--colsep', "\\t"]
    
    # Get the number of available nodes and CPUs per node
    numNodes = IrgPbsFunctions.getNumNodesInList(options.nodesListPath)
    
    # We assume all machines have the same number of CPUs (cores)
    cpusPerNode = IrgSystemFunctions.get_num_cpus()
    
    # We don't actually know the best number here!
    threadsPerCpu = 4
    
    # Set the optimal number of processes if the user did not specify
    if not options.numProcesses:
        options.numProcesses = numNodes * cpusPerNode * threadsPerCpu
        
    # Note: mapproject can run with multiple threads on non-ISIS data but we don't use that
    #       functionality here since we call mapproject with one tile at a time.
        
    # No need for more processes than their are tiles!
    if options.numProcesses > numTiles:
        options.numProcesses = numTiles
    
    # Build the command line that will be passed to GNU parallel
    # - The numbers in braces will receive the values from the text file we wrote earlier
    # - The output path used here does not matter since spawned copies compute the correct tile path.
    commandList   = ['parallel_mapproject.py',  '--pixelStartX', '{1}',
                                                     '--pixelStartY', '{2}',
                                                     '--pixelStopX',  '{3}',
                                                     '--pixelStopY',  '{4}',
                                                     '--threads', '1', # Only use on thread internally, parallel will handle things.
                                                     '--work-dir', tempFolder,
                                                     options.demPath, options.imagePath, options.outputPath]
    if options.convertTiles:
        commandList = commandList + ['--convert-tiles']
    if options.suppressOutput:
        commandList = commandList + ['--suppress-output']
    commandList   = commandList + options.extraArgs # Append other options
    commandString = IrgSystemFunctions.argListToString(commandList)    
    
    # Use GNU parallel call to distribute the work across computers
    # - This call will wait until all processes are finished
    IrgPbsFunctions.runInGnuParallel(options.numProcesses, commandString, argumentFilePath, parallelArgs, options.nodesListPath, not options.suppressOutput)


    # Build a gdal VRT file which is composed of all the processed tiles
    vrtPath = os.path.join(tempFolder, 'mosaic.vrt')
    cmd = "gdalbuildvrt -resolution highest " + vrtPath + " " + tempFolder + "*_.tif";
    print cmd
    os.system(cmd)
    
    # Convert VRT file to final output file
    cmd = "gdal_translate -co compress=lzw " + vrtPath + " " + options.outputPath;
    print cmd
    os.system(cmd)
    #IrgSystemFunctions.executeCommand(cmd, False, True)
    
    # Clean up temporary files
    if not options.keep:
        IrgFileFunctions.removeFolderIfExists(tempFolder)


    endTime = time.time()

    print "Finished in " + str(endTime - startTime) + " seconds."

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
