#!/bin/sh

--> DO NOT RUN THIS ONE!


## Script expects data to be organized by date and location
## i.e.  AN_2009_10_16  

## ROOTDIR should be set to the location where all of the
## 'data-day' directories are as shown above
## This script should also be placed in 'ROOTDIR'

SOFTWAREDIR=/u/smcmich1/icebridge/upload_software
ROOTDIR=$1

## Test on single campaign? 
## for DATADIR in `ls | egrep "GR_2015"`
for DATADIR in `ls $ROOTDIR | egrep "AN_20|GR_20"`
do

echo "Looking up details for folder: $DATADIR"

j=`echo $DATADIR |cut -f1 -d'.'`

#USED FOR 2015 GR WHERE SPRING AND FALL CAMPAIGNS WERE FLOWN
MONTH=`echo $DATADIR | cut -f2 -d'.'`
case $j in
AN_2009*)
CAMP=2009_AN_NASA
PLATFORM=DC-8
TAILNUM=N817NA
;;
AN_2010*)
CAMP=2010_AN_NASA
PLATFORM=DC-8
TAILNUM=N817NA
;;
AN_2011*)
CAMP=2011_AN_NASA
PLATFORM=DC-8
TAILNUM=N817NA
;;
AN_2012*)
CAMP=2012_AN_NASA
PLATFORM=DC-8
TAILNUM=N817NA
;;
AN_2013*)
CAMP=2013_AN_NASA
PLATFORM=P-3B
TAILNUM=N426NA
;;
AN_2014*)
CAMP=2014_AN_NASA
PLATFORM=DC-8
TAILNUM=N817NA
;;
AN_2016*)
CAMP=2016_AN_NASA
PLATFORM=DC-8
TAILNUM=N817NA
;;
GR_201[0-4]*)
CAMP=2010_GR_NASA
PLATFORM=DC-8
TAILNUM=N817NA
;;
GR_2015*)
CAMP=2015_GR_NASA
if [ "$MONTH" == "03" ] || [ "$MONTH" == "04" ] || [ "$MONTH" == "05" ]
then
  PLATFORM=C-130
  TAILNUM=N439NA
else
  PLATFORM=HU-25C
  TAILNUM=N525NA
fi
;;
GR_2016*)
CAMP=2016_GR_NASA
PLATFORM=HU-25C
TAILNUM=N525NA
;;
*)
echo "THE CAMPAIGN, PLATFORM AND TAIL NUMBER COULD NOT BEFOUND FOR $DATADIR"
;;
esac


#cd ${ROOTDIR}/${DATADIR}
echo ""
echo "---Running New Batch Directory ---"
echo "${SOFTWAREDIR}/IODEM3PremetSpatialMaker_fast.sh $CAMP $PLATFORM $TAILNUM ${ROOTDIR}/${DATADIR}"
echo ""
${SOFTWAREDIR}/IODEM3PremetSpatialMaker_fast.sh $CAMP $PLATFORM $TAILNUM ${ROOTDIR}/${DATADIR}


done


