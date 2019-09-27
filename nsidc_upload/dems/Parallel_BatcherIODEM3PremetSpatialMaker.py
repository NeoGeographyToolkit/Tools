#!/bin/python

import os, sys
import extract_oib_flight_metadata

# Script expects data to be organized by date and location
# i.e.  AN_2009_10_16  

SOFTWARE_DIR = '/u/smcmich1/icebridge/upload_software'

def main(argsIn):

    rootDir  = argsIn[0] # The top level folder where the other flights are
    yyyymmdd = argsIn[1] # The date

    # Check all subfolders
    subfolders = os.listdir(rootDir)
    for f in subfolders:

        # Skip other folders
        if ( (('AN_20' not in f) and ('GR_20' not in f)) or (yyyymmdd not in f)
             or ('cameras' in f) or ('summaries' in f)):
            print 'Skipping folder ' + f
            continue

        print 'Looking up details for folder: ' + f

        # Get the flight name
        try:
            parts = f.split('_')
            name  = parts[0] + '_' + parts[1]
        except:
            print 'Failed to parse folder: ' + f
            continue

        # Look up the flight info
        try:
            (campaign, platform, tailnum) = extract_oib_flight_metadata.main([name])
        except Exception as e:
            print str(e)
            print "THE CAMPAIGN, PLATFORM AND TAIL NUMBER COULD NOT BE FOUND FOR " + name
            continue


        # Generate files
        cmd = ('python ' + SOFTWARE_DIR + '/parallel_premet_maker.py --camp "' + campaign +
               '" --platform "' + platform + '" --tailnum "' + tailnum +
               '" --folder ' + os.path.join(rootDir, f))
        print '\n---Running New Batch Directory ---'
        print cmd +'\n'
        os.system(cmd)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
