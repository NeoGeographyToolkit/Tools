#!/usr/bin/python

#!/usr/bin/env python


# Generate a DEM for a single pair of icebridge images

import os, sys, optparse, datetime
sys.path.insert(0,  os.environ['HOME'] + '/projects/StereoPipeline/src/asp/Python')

# For sparse_disp
os.environ["ASP_PYTHON_MODULES_PATH"] = os.environ['HOME'] + '/projects/BinaryBuilder/StereoPipelinePythonModules/lib64/python2.6/site-packages:' + os.environ['HOME'] + '/projects/BinaryBuilder/StereoPipelinePythonModules/lib64/python2.6/site-packages/GDAL-1.10.0-py2.6-linux-x86_64.egg/osgeo:' + os.environ['HOME'] + '/projects/BinaryBuilder/StereoPipelinePythonModules/lib'

import asp_system_utils, asp_alg_utils, asp_geo_utils


def parseDateTimeStrings(dateString, timeString):
    '''Parse strings in the format 20110323_17433900'''
    
    MILLISECOND_TO_MICROSECOND = 10000
    
    #print dateString
    #print timeString
    
    year    = int(dateString[0:4])
    month   = int(dateString[4:6])
    day     = int(dateString[6:8])
    hour    = int(timeString[0:2])
    minute  = int(timeString[2:4])
    second  = int(timeString[4:6])
    usecond = 0
    if len(timeString) > 6:
        usecond = int(timeString[6:8]) * MILLISECOND_TO_MICROSECOND
    
    return datetime.datetime(year, month, day, hour, minute, second, usecond)

def findMatchingLidarFile(imageFile, lidarFolder):
    '''Given an image file, find the best lidar file to use for alignment.'''
    
    print imageFile
    
    # Get the image date time
    imageName       = os.path.basename(imageFile)
    parts           = imageName.replace('.tif','').split('_')
    imageDateString = parts[1]
    imageTimeString = parts[2]
    imageDateTime   = parseDateTimeStrings(imageDateString, imageTimeString)

    # Search for the matching file in the lidar folder.
    # - We are looking for the closest lidar time that starts BEFORE the image time.
    # - It is possible for an image to span lidar files, we will address that if we need to!
    bestTimeDelta = datetime.timedelta.max
    bestLidarFile = 'NA'
    lidarFiles    = os.listdir(lidarFolder)
    zeroDelta     = datetime.timedelta()
    for f in lidarFiles:
        if '.csv' not in f: # Skip other files
            continue
        # Extract time for this file
        lidarPath       = os.path.join(lidarFolder, f)
        parts           = f.replace('.','_').split('_')
        lidarDateString = parts[1]
        lidarTimeString = parts[2]
        lidarDateTime   = parseDateTimeStrings(lidarDateString, lidarTimeString)
        
        # Compare time to the image time
        timeDelta       = imageDateTime - lidarDateTime
        if ( (timeDelta > zeroDelta) and (timeDelta < bestTimeDelta) ):
            bestLidarFile = lidarPath
            bestTimeDelta = timeDelta

    #print 'Found matching lidar file ' + bestLidarFile
    #print bestTimeDelta

    if bestLidarFile == 'NA':
        raise Exception('Failed to find matching lidar file for image ' + imageFile)

    return bestLidarFile


def main(argsIn):

    try:
        usage = 'usage: process_icebridge_pair.py <imageA> <imageB> <cameraA> <cameraB> <lidar_folder> <output_folder>[--help]\n  '
        parser = optparse.OptionParser(usage=usage)

        parser.add_option('--south', action='store_true', default=False, dest='isSouth',  
                          help='MUST be set if the images are in the southern hemisphere.')
                          
        parser.add_option('--lidar-overlay', action='store_true', default=False, dest='lidarOverlay',  
                          help='Generate a lidar overlay for debugging.')

        parser.add_option('--bundle-adjust', action='store_true', default=False, dest='bundleAdjust',  
                          help='Run bundle adjustment between the two images.')

        parser.add_option('--num-threads', dest='numThreads', default=None,
                          type='int', help='The number of threads to use for processing.')

        parser.add_option('--dem-resolution', dest='demResolution', default=0.4,
                          type='float', help='Generate output DEMs at this resolution.')

        parser.add_option('--align-max-displacement', dest='maxDisplacement', default=20,
                          type='float', help='Max displacement value passed to pc_align.')

        parser.add_option('--use-sgm', action='store_true', default=False, dest='use_sgm',  
                          help='If to use SGM.')

        parser.add_option('--pc-align', action='store_true', default=False, dest='pc_align',  
                          help='If to use pc_align.')

        (options, args) = parser.parse_args(argsIn)

        if len(args) < 7:
            print usage
            return 0

        imageA       = args[1]
        imageB       = args[2]
        cameraA      = args[3]
        cameraB      = args[4]
        lidarFolder  = args[5]
        outputFolder = args[6]

    except optparse.OptionError, msg:
        raise Usage(msg)

    # Pick the output projection to be used
    PROJ_STRING_NORTH = '"+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"'
    PROJ_STRING_SOUTH = '"+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"'
    
    projString = PROJ_STRING_NORTH
    if options.isSouth:
        projString = PROJ_STRING_SOUTH

    # Check the inputs
    for f in [imageA, imageB, cameraA, cameraB, lidarFolder]:   
        if not os.path.exists(f):
            print 'Input file '+ f +' does not exist!'
            return 0
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)

    print 'Searching for matching lidar file...'
    lidarFile = findMatchingLidarFile(imageA, lidarFolder)
    print 'Found matching lidar file ' + lidarFile

    # Does this ever change?
    csvFormatString = '"1:lat 2:lon 3:height_above_datum"'

    suppressOutput = False
    redo           = False

    print '\nStarting processing...'

    outputPrefix  = os.path.join(outputFolder, 'out')
        
    baseArgString = ('%s %s %s %s %s -t nadirpinhole --alignment-method epipolar' 
                     % (imageA, imageB, cameraA, cameraB, outputPrefix))
    
    threadText = ''
    if options.numThreads:
        threadText = ' --threads ' + str(options.numThreads) +' '
    
    # BUNDLE_ADJUST
    if options.bundleAdjust:
        # TODO: Solve for intrinsics?
        bundlePrefix = os.path.join(outputFolder, 'bundle/out')
        cmd = ('bundle_adjust %s %s %s %s -o %s %s -t nadirpinhole --local-pinhole' 
                     % (imageA, imageB, cameraA, cameraB, bundlePrefix, threadText))
        print(cmd)
        # Point to the new camera models
        cameraA = bundlePrefix +'-'+ os.path.basename(cameraA)
        cameraB = bundlePrefix +'-'+ os.path.basename(cameraB)
        asp_system_utils.executeCommand(cmd, cameraA, suppressOutput, redo)

        # Update the baseArgString
        baseArgString = ('%s %s %s %s %s -t nadirpinhole --alignment-method epipolar' 
                         % (imageA, imageB, cameraA, cameraB, outputPrefix))
    
    if options.use_sgm:
        # PPRC
        cmd = ('stereo_pprc %s %s' % (baseArgString, threadText))
        pprcOutput = outputPrefix + '-L.tif'
        asp_system_utils.executeCommand(cmd, pprcOutput, suppressOutput, redo)

        # CORR
        # - This should be single threaded to use the SGM processing.
        # - TODO: Use blob filtering to reduce outliers?
        # - TODO: Can we shrink the search range?
        correlationArgString = ('--threads 1 --xcorr-threshold -1 --corr-kernel 7 7' 
                               +' --corr-tile-size 6400 --cost-mode 4 --use-sgm --subpixel-mode 0')
                               #+ ' --corr-blob-filter 100 --compute-low-res-disparity-only')
        cmd = ('stereo_corr %s %s ' % (correlationArgString, baseArgString))
        corrOutput = outputPrefix + '-D.tif'
        asp_system_utils.executeCommand(cmd, corrOutput, suppressOutput, redo)

        #raise Exception('DEBUG')

        # RFNE
        cmd = ('stereo_rfne --subpixel-mode 0 %s %s' % (baseArgString, threadText))
        rfneOutput = outputPrefix + '-RD.tif'
        asp_system_utils.executeCommand(cmd, rfneOutput, suppressOutput, redo)

        # FLTR
        filterArgString = '--rm-cleanup-passes 0 --median-filter-size 5 ' + \
                          '--texture-smooth-size 17 --texture-smooth-scale 0.14'
        cmd = ('stereo_fltr %s %s %s' % (filterArgString, baseArgString, threadText))
        fltrOutput = outputPrefix + '-F.tif'
        asp_system_utils.executeCommand(cmd, fltrOutput, suppressOutput, redo)

        # TRI
        cmd = ('stereo_tri %s %s' % (baseArgString, threadText))
        triOutput = outputPrefix + '-PC.tif'
        asp_system_utils.executeCommand(cmd, triOutput, suppressOutput, redo)

        #raise Exception('DEBUG')
    else: # No SGM
        cmd = ('stereo %s %s --corr-blob-filter 10000 --subpixel-mode 3 --corr-timeout 120' % \
               (baseArgString, threadText))
        
        # Fine level control when using subpixel-mode 3. TODO: This
        # needs to be done better.
        #lOutput = outputPrefix + '-L.tif'
        #asp_system_utils.executeCommand(cmd, lOutput, suppressOutput, redo)

        #cmd = ('stereo -e 1 --stop-point 2 --compute-low-res-disparity-only %s %s --corr-max-levels 2 --corr-seed-mode 3 --subpixel-mode 3' % (baseArgString, threadText)) 
        #dOutput = outputPrefix + '-D_sub.tif'
        #asp_system_utils.executeCommand(cmd, dOutput, suppressOutput, redo)

        #cmd = ('stereo_corr --skip-low-res-disparity-comp --corr-blob-filter 0 %s %s --corr-max-levels 2 --corr-seed-mode 3 --subpixel-mode 3' % (baseArgString, threadText))
        #dOutput = outputPrefix + '-D.tif'
        #asp_system_utils.executeCommand(cmd, dOutput, suppressOutput, redo)

        #cmd = ('stereo -e 2 --skip-low-res-disparity-comp --corr-blob-filter 0 %s %s --corr-max-levels 2 --corr-seed-mode 3 --subpixel-mode 3' % (baseArgString, threadText))
        
        triOutput = outputPrefix + '-PC.tif'
        asp_system_utils.executeCommand(cmd, triOutput, suppressOutput, redo)

    # point2dem on the result of ASP
    cmd = ('point2dem --tr %lf --t_srs %s %s %s --errorimage' 
           % (options.demResolution, projString, triOutput, threadText))
    p2dOutput = outputPrefix + '-DEM.tif'
    asp_system_utils.executeCommand(cmd, p2dOutput, suppressOutput, redo)

    if options.pc_align:
        # PC_ALIGN
        alignPrefix = os.path.join(outputFolder, 'align/out')
        alignOptions = ( ('--max-displacement %f --csv-format %s ' +   \
                          '--save-inv-transformed-reference-points') % \
                         (options.maxDisplacement, csvFormatString))
        cmd = ('pc_align %s %s %s -o %s %s' %
               (alignOptions, triOutput, lidarFile, alignPrefix, threadText))
        alignOutput = alignPrefix+'-trans_reference.tif'
        asp_system_utils.executeCommand(cmd, alignOutput, suppressOutput, redo)
        
        # POINT2DEM on the aligned PC file
        cmd = ('point2dem --tr %lf --t_srs %s %s %s --errorimage' 
               % (options.demResolution, projString, alignOutput, threadText))
        p2dOutput = alignPrefix+'-trans_reference-DEM.tif'
        asp_system_utils.executeCommand(cmd, p2dOutput, suppressOutput, redo)
       
    # Create a symlink to the DEM in the main directory
    demSymlinkPath = os.path.join(outputFolder, 'DEM.tif')
    print("ln -s " + os.path.abspath(p2dOutput) + " " + demSymlinkPath)
    os.symlink(os.path.abspath(p2dOutput), demSymlinkPath)

    cmd = ('geodiff --absolute --csv-format %s %s %s -o %s' % \
           (csvFormatString, p2dOutput, lidarFile, outputPrefix))
    print(cmd)
    asp_system_utils.executeCommand(cmd, outputPrefix + "-diff.csv", suppressOutput, redo)

    # HILLSHADE
    hillOutput = outputPrefix+'-DEM_HILLSHADE.tif'
    cmd = 'hillshade ' + p2dOutput +' -o ' + hillOutput
    asp_system_utils.executeCommand(cmd, hillOutput, suppressOutput, redo)
    
    # COLORMAP
    colormapMin = -10
    colormapMax =  10
    colorOutput = outputPrefix+'-DEM_CMAP.tif'
    cmd = ('colormap --min %f --max %f %s -o %s' 
           % (colormapMin, colormapMax, p2dOutput, colorOutput))
    asp_system_utils.executeCommand(cmd, colorOutput, suppressOutput, redo)

    if options.lidarOverlay:
        LIDAR_DEM_RESOLUTION     = 5
        LIDAR_PROJ_BUFFER_METERS = 100
    
        # Get buffered projection bounds of this image
        demGeoInfo = asp_geo_utils.getImageGeoInfo(p2dOutput, getStats=False)
        projBounds = demGeoInfo['projection_bounds']
        minX = projBounds[0] - LIDAR_PROJ_BUFFER_METERS # Expand the bounds a bit
        minY = projBounds[2] - LIDAR_PROJ_BUFFER_METERS
        maxX = projBounds[1] + LIDAR_PROJ_BUFFER_METERS
        maxY = projBounds[3] + LIDAR_PROJ_BUFFER_METERS

        # Generate a DEM from the lidar point cloud in this region        
        lidarDemPrefix = os.path.join(outputFolder, 'cropped_lidar')
        cmd = ('point2dem --t_projwin %f %f %f %f --tr %lf --t_srs %s %s %s --csv-format %s -o %s' 
               % (minX, minY, maxX, maxY,
                  LIDAR_DEM_RESOLUTION, projString, lidarFile, threadText, 
                  csvFormatString, lidarDemPrefix))
        lidarDemOutput = lidarDemPrefix+'-DEM.tif'
        asp_system_utils.executeCommand(cmd, lidarDemOutput, suppressOutput, redo)
            
        colorOutput = lidarDemPrefix+'-DEM_CMAP.tif'
        cmd = ('colormap --min %f --max %f %s -o %s' 
               % (colormapMin, colormapMax, lidarDemOutput, colorOutput))
        asp_system_utils.executeCommand(cmd, colorOutput, suppressOutput, redo)

    print 'Finished!'

# Run main function if file used from shell
if __name__ == "__main__":
    sys.exit(main(sys.argv))



