#!/bin/sh

#####
# MJS 2018 03 05
#
# This script should be ran from the directory full of IODEM3 tif files to parse are.
# (Or the DIR variable can be set, instead of DIR=pwd)
# The PLATFORM, TAILNUM and CAMPAIGN vars in the section need to be set for each campaign.
#
# Data files must start IODEM3 and end with .tif to be found and parsed.
##
echo ""
echo "Did you Set Campaign, Tail Number and Platform in this scirpt or on CLI?"
echo "Usage: IODEM3PremetSpatialMaker.sh <CAMPAIGN> <PLATFORM> <TAILNUMBER> <DATADIR>"

#sleep 2
#
## When using BATCHER GET CAMP PLAT AND TAILN FROM CLI 
CAMP=$1
PLATFORM=$2
TAILNUM=$3
DIR=$4

# DIR IS SET TO DATA DIRECTORY TO PROCESS
#DIR=`pwd`

##### Platfor Tail Num options
#PLATFORM=P-3B
#TAILNUM=N426NA
#PLATFORM=DC-8
#TAILNUM=N817NA
#PLATFORM=HU-25C
#TAILNUM=N525NA
#PLATFORM=C-130
#TAILNUM=N439NA
#PLATFORM=HU-25C
#TAILNUM=N525NA

GDALDIR=/u/smcmich1/programs/StereoPipeline-2.6.0-2018-02-06-x86_64-Linux/bin

### END CONFIG

create_spo (){
UL_IS_W=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Upper Left" | cut -f3 -d'(' | cut -f1 -d','| cut -f2 -d'"'`
UL_IS_S=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Upper Left" | cut -f3 -d'(' | cut -f2 -d','| cut -f2 -d'"' | cut -f1 -d')'`
UL_LON=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Upper Left" | cut -f3 -d'(' | cut -f1 -d','| sed 's: :0:g'| awk '{degs=substr($1,1,3);mins=substr($1,5,2)/60;secs=substr($1,8,5)/3600; outp=degs+mins+secs; print(outp)}'`
UL_LAT=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Upper Left" | cut -f3 -d'(' | cut -f2 -d','| sed 's: :0:g'| awk '{degs=substr($1,1,3);mins=substr($1,5,2)/60;secs=substr($1,8,5)/3600; outp=degs+mins+secs; print(outp)}'`
if [ "$UL_IS_W" == "W" ]
 then
 UL_LON=-$UL_LON
fi
if [ "$UL_IS_S" == "S" ]
 then
 UL_LAT=-$UL_LAT
fi

UR_IS_W=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Upper Right" | cut -f3 -d'(' | cut -f1 -d','| cut -f2 -d'"'`
UR_IS_S=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Upper Right" | cut -f3 -d'(' | cut -f2 -d','| cut -f2 -d'"' | cut -f1 -d')'`
UR_LON=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Upper Right" | cut -f3 -d'(' | cut -f1 -d','| sed 's: :0:g'| awk '{degs=substr($1,1,3);mins=substr($1,5,2)/60;secs=substr($1,8,5)/3600; outp=degs+mins+secs; print(outp)}'`
UR_LAT=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Upper Right" | cut -f3 -d'(' | cut -f2 -d','| sed 's: :0:g'| awk '{degs=substr($1,1,3);mins=substr($1,5,2)/60;secs=substr($1,8,5)/3600; outp=degs+mins+secs; print(outp)}'`
if [ "$UR_IS_W" == "W" ]
 then
 UR_LON=-$UR_LON
fi
if [ "$UR_IS_S" == "S" ]
 then
 UR_LAT=-$UR_LAT
fi

LR_IS_W=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Lower Right" | cut -f3 -d'(' | cut -f1 -d','| cut -f2 -d'"'`
LR_IS_S=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Lower Right" | cut -f3 -d'(' | cut -f2 -d','| cut -f2 -d'"' | cut -f1 -d')'`
LR_LON=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Lower Right" | cut -f3 -d'(' | cut -f1 -d','| sed 's: :0:g'| awk '{degs=substr($1,1,3);mins=substr($1,5,2)/60;secs=substr($1,8,5)/3600; outp=degs+mins+secs; print(outp)}'`
LR_LAT=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Lower Right" | cut -f3 -d'(' | cut -f2 -d','| sed 's: :0:g'| awk '{degs=substr($1,1,3);mins=substr($1,5,2)/60;secs=substr($1,8,5)/3600; outp=degs+mins+secs; print(outp)}'`
if [ "$LR_IS_W" == "W" ]
 then
 LR_LON=-$LR_LON
fi
if [ "$LR_IS_S" == "S" ]
 then
 LR_LAT=-$LR_LAT
fi

LL_IS_W=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Lower Left" | cut -f3 -d'(' | cut -f1 -d','| cut -f2 -d'"'`
LL_IS_S=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Lower Left" | cut -f3 -d'(' | cut -f2 -d','| cut -f2 -d'"' | cut -f1 -d')'`
LL_LON=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Lower Left" | cut -f3 -d'(' | cut -f1 -d','| sed 's: :0:g'| awk '{degs=substr($1,1,3);mins=substr($1,5,2)/60;secs=substr($1,8,5)/3600; outp=degs+mins+secs; print(outp)}'`
LL_LAT=`${GDALDIR}/gdalinfo ${DIR}/${FILENAME}  | grep "Lower Left" | cut -f3 -d'(' | cut -f2 -d','| sed 's: :0:g'| awk '{degs=substr($1,1,3);mins=substr($1,5,2)/60;secs=substr($1,8,5)/3600; outp=degs+mins+secs; print(outp)}'`
if [ "$LL_IS_W" == "W" ]
 then
 LL_LON=-$LL_LON
fi
if [ "$LL_IS_S" == "S" ]
 then
 LL_LAT=-$LL_LAT
fi

echo $UL_LON $UL_LAT > ${DIR}/${FILENAME}.spo
echo $UR_LON $UR_LAT >> ${DIR}/${FILENAME}.spo
echo $LR_LON $LR_LAT >> ${DIR}/${FILENAME}.spo
echo $LL_LON $LL_LAT >> ${DIR}/${FILENAME}.spo

}


create_premet (){
  DATE=`echo ${FILENAME} | cut -f2 -d'_' | awk '{year=substr($1,1,4);mon=substr($1,5,2);day=substr($1,7,2); print(year,"-",mon,"-",day)}' | sed $1 's: ::g'`
  TIME=`echo ${FILENAME} | cut -f3 -d'_' | awk '{hr=substr($1,1,2);min=substr($1,3,2);sec=substr($1,5,2);decsec=substr($1,7,4); print(hr,":",min,":",sec,".",decsec)}' | sed $1 's: ::g'`
  echo "VersionID_local=001" > ${DIR}/${FILENAME}.premet
  echo "Begin_date=$DATE" >> ${DIR}/${FILENAME}.premet
  echo "End_date=$DATE" >> ${DIR}/${FILENAME}.premet
  echo "Begin_time=$TIME" >> ${DIR}/${FILENAME}.premet
  echo "End_time=$TIME" >> ${DIR}/${FILENAME}.premet
  echo "CampaignShortName=$CAMP" >> ${DIR}/${FILENAME}.premet
  echo "Container=AdditionalAttributes" >> ${DIR}/${FILENAME}.premet
  echo "AdditionalAttributeName=ThemeID" >> ${DIR}/${FILENAME}.premet
  echo "ParameterValue=$CAMP" >> ${DIR}/${FILENAME}.premet
  echo "Container=AdditionalAttributes" >> ${DIR}/${FILENAME}.premet
  echo "AdditionalAttributeName=AircraftID" >> ${DIR}/${FILENAME}.premet
  echo "ParameterValue=$TAILNUM" >> ${DIR}/${FILENAME}.premet
  echo "Container=AssociatedPlatformInstrumentSensor" >> ${DIR}/${FILENAME}.premet
  echo "AssociatedPlatformShortName=$PLATFORM" >> ${DIR}/${FILENAME}.premet
  echo "AssociatedInstrumentShortName=DMS" >> ${DIR}/${FILENAME}.premet
  echo "AssociatedSensorShortName=DMS" >> ${DIR}/${FILENAME}.premet
}



echo "Writing Files to $DIR"
#for FILENAME in `ls $DIR | grep "^IODEM3" | grep ".tif$"`
for FILENAME in `ls $DIR  | grep ".tif$"`
do
  echo IODEM3PremetSpatialMaker.sh is running ${DIR}/${FILENAME}
  create_spo
  create_premet

done

