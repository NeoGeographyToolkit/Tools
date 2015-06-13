
import sys
import os


'''
    Script to copy the metadata from a GeoTiff file to a plain Tiff file
'''

def copyGeoTiffInfo(inputGeoTiffPath, inputTiffPath, outputPath):
    '''Copies geo information from a GeoTiff to a plain Tiff'''

    # Extract geo information from the input geotiff
    geoInfoPath = inputGeoTiffPath + '.metadata'
    cmd = 'listgeo ' + inputGeoTiffPath +' > '+ geoInfoPath
    print cmd
    os.system(cmd)
    
    # Paste the geo information into the output file
    cmd = 'geotifcp ' + inputTiffPath +' '+ outputPath +' -g '+ geoInfoPath
    print cmd
    os.system(cmd)
    
    # Clean up temporary file
    os.remove(geoInfoPath)

def main():
    
    # Get the input arguments from the command line
    if len(sys.argv) < 3:
        raise Exception('Usage: copyGeotiffInfo.py <input geotiff path> <input tiff path> [output geotiff path]')
    
    inputGeoTiffPath = sys.argv[1]
    inputTiffPath    = sys.argv[2]
    if len(sys.argv) >= 4:
        outputPath = sys.argv[3]
    else:
        outputPath = inputGeoTiffPath + '_geo.tif'
    
    # Call function to do the work
    copyGeoTiffInfo(inputGeoTiffPath, inputTiffPath, outputPath)    
    
    
    print 'Finished copying metadata to output file: ' + outputPath


if __name__ == "__main__":
    sys.exit(main())