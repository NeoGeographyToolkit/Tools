 You will need jave 1.8.54 or newer, and gdalinfo to use this.

 If you download that and unpackage it you'll see:

1) some directories named AN_<DATE> and GR_<DATE>
which contain 3 sample files from each of those dates
so you can see how the scripts work.


2) IODEM3PremetSpatialMaker.sh script. Which generates .premet and .spo
files back into each data dir. These are used by our SIPSMetGen tool.

This script takes flags, or can be edited to hard code
these values. Run as below unless edited:

./IODEM3PremetSpatialMaker.sh <CAMPAIGN> <PLATFORM> <TAILNUMBER> <DATA_DIRECTORY_TO_PROCESS>

It is commented well enough I think, and has samples of what those values
should be as comments in it.

Last, in this script
GDALDIR  variable must be set to the the path to the 'gdalinfo' binary.
(i.e. GDALDIR=/usr/bin/  , if gdalinfo is at /usr/bin/gdalinfo)



3) A BatcherIODEM3PremetSpatialMaker.sh scirpt.
This script will run the IODEM3PremetSpatialMaker.sh script above in a
batch on all of the directories it finds that start with "AN_ and/or GR_"


It is also commented well enough.  It will determine the campaign,
then lookup the correct platform and tail number for each campaign,
and call the above script for the correct data dir for each data dir
it finds. (Has values for 2009_AN_NASA - 2016_AN_NASA and 2010_GR_NASA
-2016_GR_NASA, feel free to add new ones if you have more campaigns,
or let me know which ones to add and I can update it.)

4) MetGen.Start Script, which will allow you to run
SIPSMetGen (or SMG) for a directory of data, assuming
the paths in the Primary.cfg file are set.

 MetGen.Start has a JAVA_HOME variable you will need
to set to where the java binary (1.8.54 or newer)
is on your system.

Run it with:  MetGen.Start -config ./Primary.cfg

5) Primary.cfg file, which has paths that you will need
to edit.

  If you replace the each of the paths from '=' to '/ScriptTest/'
with the location you unpackage the .tgz file which contains the
ScriptTest dir you should be good to go using the batch script
described below.  You'll need to change 'DATADIR' as well if you  
run it one at a time. I'd recommend changing the batch script
to have a list of only 1 directory instead, and let it run
the single dir. This will make sense when you get to #7 below.

6) IODEM3#001.MCF, valids.cfg and push.cfg which you won't really need to
do anything to but are needed by SMG.

7) SMG_IODEM3Batcher.sh  Script runs SMG on all the  
directories listed in the for loop. It replaces 'DATADIR'
in the Primary config file, with the correct dir from
the for loop list, and runs through the whole list  
and sets the DATADIR back to be run again.


---
I know this sounds like a lot, but it's not bad once you dig into it.

Essentially you need to Edit the IDOEM3PremetSpatialMaker.sh script
on each data dir. Manually or using the Batch script.

Then, run MetGen on each data directory, Manually or using
it's batch script.
---
PremetSpatialMaker
---
Make sure you've set GDALDIR in the IODEM3PremetSpatialMaker.sh
and are feeding it good values through the CLI or from the
Batch script for campaign, platform, tailnumber, and  data_directory.

Make sure the ROOTDIR is set properly if you are using the
BatcherIODEM3PremetSpatialMaker.sh script.
---
SMG
---
Edit the Primary.cfg File before running SMG. Set paths as
noted above.

Edit the MetGen.Start script to set JAVA_HOME

Then run it either Manually or using the SMG_IODEM3Batcher.sh.

For the SMG_IODEM3Batcher.sh, set ROOTDIR, and set
the list of directories to process.  The one I have
has a 'for i in ...' commented out with a bunch
of GR_<DATE> and AN_<DATE> values, and a short
one with only two data dirs to test it with.
---

Feel free to ask any questions. I'm happy to help whenever I can.

schwabm@nsidc.org
