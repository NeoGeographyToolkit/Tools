#!/bin/sh

KMLFILE=doc.kml

#Find  all RDSISC(O)4 tif files in current directory
for FILE in `find ./ -maxdepth 1 -name "RDSISC*tif"`
do
  #Get  segment to match  in KML file
  SEGMENT=`echo $FILE  | cut -f2-5 -d'_'`

  LINES=`egrep -A4 $SEGMENT ./$KMLFILE`
  RES=$?
  if [ $RES -eq "0" ]
  then

      #FOR EACH MATCH SPLIT TEXT  and  Grap  LON  LATS
     LINES=`egrep -A4 $SEGMENT ./$KMLFILE`
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
        echo  "$BGLAT $SMLON" >  ${FILE}.spo
        echo  "$BGLAT $BGLON" >>  ${FILE}.spo
        echo  "$SMLAT $BGLON" >>  ${FILE}.spo
        echo  "$SMLAT $SMLON" >>  ${FILE}.spo
     else
        echo  "$BGLAT $BGLON" >  ${FILE}.spo
        echo  "$BGLAT $SMLON" >>  ${FILE}.spo
        echo  "$SMLAT $SMLON" >>  ${FILE}.spo
        echo  "$SMLAT $BGLON" >>  ${FILE}.spo

     fi
  else
    echo "ERROR *****  $FILE ****  NOT  FOUND IN  KML FILE  $KMLFILE ****  "
  fi

done


