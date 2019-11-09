#!/bin/sh


## FINDS ALL RDSISC4 DATA FILES AND RUNS RDSISC4_Parser.py on them
##
## Copy this script and the parser script into the directory to create premet and spo files for
## 

PREFIX=RDSISC4_
SUFFIX=_ortho.tif
SCRIPT=RDSISC4_Parser.py

PWD=`pwd`

for FILE in `ls $PWD | grep ^$PREFIX | grep ${SUFFIX}$`
do
    #echo "${PWD}/${SCRIPT} $FILE"
    ${PWD}/${SCRIPT} $FILE
done



