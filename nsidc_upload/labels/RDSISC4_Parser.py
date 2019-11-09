#!/tools/python/bin/python

import os
import numpy as np
import subprocess
import sys


## SCRIPT TAKES FILE NAME IS ARGUMENT. FAILS IF FILE NOT FOUND
## Uses GDALINFO WHICH MUST BE INSTALLED ON SYSTEM TO RUN and BE IN YOUR PATH
## PARSES CORNER POINTS PUT THESE IN .SPO FILE FOR GIVEN FILE
##
## CREATES PREMET FILE FOR WHOLE DAY AS TIME RANGE FOR DATE IN FILE NAME
## FILE NAME FORMAT SHOULD BE: RDSISC4_2017_04_11_00008_classified_ortho.tif
## PARTICULARLY THE YEAR MONTH AND DAY MUST BE IN PLACES 2 3 AND 4 SPLIT BY _
##
## SCRIPT USES INPUT VALUES SET BELOW TO INSERT AIRCRAFT, CAMPAIGN AND INSTRUMENT INFO
## MUST SET THESE FOR EACH OIB CAMPAIGN!!!!!!!!

#### SET THESE FOR EACH CAMPAIGN CORRECTLY !!!!
instrument='DMS'
sensor='DMS'

#local version id probably it is 001 
lvid='001'
#
# just set time to whole day hopefully users won't search for data by the second
btimestamp='00:00:00'
etimestamp='23:59:59'


# GRABS LINS FROM GDALINFO AND PUTS THE COORDINATES IN LAT LON ARRAYS
def GetCornerCoordinates(FileName):
    #print (FileName)
    GdalInfo = subprocess.check_output('gdalinfo {}'.format(FileName), shell=True)
    GdalInfo = GdalInfo.split('\n') # Creates a line by line list.
    CornerLats, CornerLons = np.zeros(4), np.zeros(4) 
    for line in GdalInfo:
        #print "LINE"
        #print line
        if line[:10] == 'Upper Left':
            CornerLats[0], CornerLons[0] = GetLatLon(line)
        if line[:10] == 'Lower Left':
            CornerLats[1], CornerLons[1] = GetLatLon(line)
        if line[:11] == 'Lower Right':
            CornerLats[2], CornerLons[2] = GetLatLon(line)
        if line[:11] == 'Upper Right':
            CornerLats[3], CornerLons[3] = GetLatLon(line)
    return CornerLats, CornerLons 

# CONVERTS TO DECIMAL DEGREES
def GetLatLon(line):
    coords = line.split(') (')[1]
    coords = coords[:-1]
    LonStr, LatStr = coords.split(',')
    # Longitude
    LonStr = LonStr.split('d')    # Get the degrees, and the rest
    LonD = int(LonStr[0])
    LonStr = LonStr[1].split('\'')# Get the arc-m, and the rest
    LonM = int(LonStr[0])
    LonStr = LonStr[1].split('"') # Get the arc-s, and the rest
    LonS = float(LonStr[0])
    Lon = LonD + LonM/60. + LonS/3600.
    if LonStr[1] in ['W', 'w']:
        Lon = -1*Lon
    # Same for Latitude
    LatStr = LatStr.split('d')
    LatD = int(LatStr[0])
    LatStr = LatStr[1].split('\'')
    LatM = int(LatStr[0])
    LatStr = LatStr[1].split('"')
    LatS = float(LatStr[0])
    Lat = LatD + LatM/60. + LatS/3600.
    if LatStr[1] in ['S', 's']:
        Lat = -1*Lat
    #print(Lat)
    return Lat, Lon


# USAGE STATEMENT
def usage(leng):
    if len(leng) < 2 :
       print "Usage : RDSISC4_Parser.py <InputFileName>"
       sys.exit()
    
    # OPEN FILE AN VOMIT IF IT DOESNT EXIST
def openfile(inputfile):
    try:
        fh = open(inputfile)
    except IOError:
        print('Error opening file: %s does not exist' % inputfile)
        sys.exit()

# GET COORDS FROM FUNCTIONS ABOVE 
# THEN WRITE SPO FILE
def writespatial():
    CornerLats, CornerLons = GetCornerCoordinates(inputfile)
    sname=inputfile+".spo"
    spatial = open(sname, 'w')
    #spatial.write('%s %s\n' % (CornerLons[0],CornerLats[0]))
    #spatial.write('%s %s\n' % (CornerLons[1],CornerLats[1]))
    #spatial.write('%s %s\n' % (CornerLons[2],CornerLats[2]))
    spatial.write('%s %s\n' % (CornerLons[3],CornerLats[3]))
    spatial.write('%s %s\n' % (CornerLons[2],CornerLats[2]))
    spatial.write('%s %s\n' % (CornerLons[1],CornerLats[1]))
    spatial.write('%s %s\n' % (CornerLons[0],CornerLats[0]))

# WRITE PREMET FILE USE VARIABLES SET AT TOP FOR FIXED INPUT
# THE TIME STAMP IS REALLY CREATION DATE  I HARD CODED THE LONG RANGE
def parsedate():
    fname = os.path.basename(inputfile)
    inputbits=fname.split('_')
    datestamp=inputbits[1]+'-'+inputbits[2]+'-'+inputbits[3]
    return datestamp
    
# WRITE PREMET FILE USE VARIABLES SET AT TOP FOR FIXED INPUT
def writepremet(campaign, platform, aircraftid):
    pname=inputfile+".premet"
    premet = open(pname, 'w')
    premet.write('VersionID_local=%s\n' % lvid)
    premet.write('Begin_date=%s\n' % datestamp)
    premet.write('Begin_time=%s\n' % btimestamp)
    premet.write('End_date=%s\n' % datestamp)
    premet.write('End_time=%s\n' % etimestamp)
    
    premet.write('Container=AssociatedPlatformInstrumentSensor\n')
    premet.write('AssociatedPlaformShortName=%s\n' % (platform))
    premet.write('AssociatedInstrumentShortName=%s\n' % (instrument))
    premet.write('AssociatedSensorShortName=%s\n' % (sensor))
    
    premet.write('Container=AdditionalAttributes\n')
    premet.write('AdditionalAttributeName=ThemeID\n')
    premet.write('ParameterValue=%s\n' % (campaign))
    
    premet.write('Container=AdditionalAttributes\n')
    premet.write('AdditionalAttributeName=AircraftID\n')
    premet.write('ParameterValue=%s\n' % (aircraftid))

#### RUN 

usage(sys.argv)
inputfile = sys.argv[1]
campaign=sys.argv[2]
platform=sys.argv[3]
aircraftid=sys.argv[4]
openfile(inputfile)
datestamp=parsedate()
writespatial()
writepremet(campaign, platform, aircraftid)


