#!/bin/sh

# Parse arguments
FILE=$1
NAME=$2
KMLFILE=$3

#Get  segment to match  in KML file
SEGMENT=`echo $NAME  | cut -f2-5 -d'_'`

LINES=`egrep -A4 $SEGMENT $KMLFILE`
RES=$?
if [ $RES -eq "0" ]
then

    #FOR EACH MATCH SPLIT TEXT  and  Grap  LON  LATS
    #LINES=`egrep -A4 $SEGMENT ./$KMLFILE`
    for LINE in `echo $LINES`
      do
        #echo  line $LINE
        TMP=`echo  $LINE | grep "latitude"`
        RESULT=$?
        if [ $RESULT -eq "0" ]
        then
        LAT=`echo $LINE | cut -f2  -d'>' |cut -f1 -d'<' `
        #else
          #echo  "ERROR ********* $FILE ********  LATITUDE NOT FOUND"
        fi
        TMP=`echo  $LINE | grep "longitude"`
        RESULT=$?
        if [ $RESULT -eq "0" ]
        then
        LON=`echo $LINE | cut -f2  -d'>' |cut -f1 -d'<' `
        fi
      done
    echo "writing ${FILE}.spo"
    echo $LON $LAT  >  ${FILE}.spo
    BGLAT=$(echo "$LAT + 0.00002" | bc)
    SMLAT=$(echo "$LAT - 0.00002" | bc)
    BGLON=$(echo "$LON + 0.00002" | bc)
    SMLON=$(echo "$LON - 0.00002" | bc)
#     echo $LAT $BGLAT

    ##  SINCE ALL IS  GR DATA JUST CHECK LONGITIUDE
    if [ $BGLON > $SMLON ]
    then
        echo  "$SMLON $BGLAT" >  ${FILE}.spo
        echo  "$BGLON $BGLAT" >>  ${FILE}.spo
        echo  "$BGLON $SMLAT" >>  ${FILE}.spo
        echo  "$SMLON $SMLAT" >>  ${FILE}.spo
     else
        echo  "$BGLON $BGLAT" >  ${FILE}.spo
        echo  "$SMLON $BGLAT" >>  ${FILE}.spo
        echo  "$SMLON $SMLAT" >>  ${FILE}.spo
        echo  "$BGLON $SMLAT" >>  ${FILE}.spo
    fi
else
  echo "ERROR *****  $FILE ****  NOT  FOUND IN  KML FILE  $KMLFILE ****  "
fi



