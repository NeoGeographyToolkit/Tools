#!/usr/bin/env python
# __BEGIN_LICENSE__
#  Copyright (c) 2009-2013, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NGT platform is licensed under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance with the
#  License. You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# __END_LICENSE__



import os, sys

SOURCE_FILE = '/u/smcmich1/icebridge/upload_software/oib_platform-aircraftId_date_summary.csv'


def main(argsIn):

    if not os.path.exists(SOURCE_FILE):
        raise Exception('Flight metadata file not found: ' + SOURCE_FILE)

    if len(argsIn) != 1:
        raise Exception('No input folder passed in!')

    folderName = argsIn[0]

    # TODO: Do we need to replace AN with AK??
    POSSIBLE_LOCATIONS = ['GR', 'AN', 'AK']
    CAMPAIGN_INDEX = 0
    PLATFORM_INDEX = 2
    AID_INDEX      = 3


    # Get information for this flight
    parts    = folderName.split('_')
    location = parts[0]
    date     = parts[1]
    #year     = date[0:4]
    #month    = date[4:6]
    #day      = date[6:8]

    if location not in POSSIBLE_LOCATIONS:
        raise Exception('Unrecognized location for folder: ' + folderName)

    match = location + '_' + date

    # Find the matching line in the data file
    matchLine = ''
    with open(SOURCE_FILE, 'r') as f:
        for line in f:
            if match in line:
                matchLine = line
                break

    if not matchLine:
        raise Exception('Folder not located in the data file!')

    # Get the data from the matching line
    parts = matchLine.split(',')
    campaign   = parts[CAMPAIGN_INDEX]
    platform   = parts[PLATFORM_INDEX]
    aircraftId = parts[AID_INDEX     ]

    return (campaign, platform, aircraftId)

if __name__ == "__main__":

    (campaign, platform, aircraftId) = main(sys.argv[1:])

    # Print out the matching data
    print 'MATCH '+ campaign +' '+ platform +' '+ aircraftId

    sys.exit(0)
