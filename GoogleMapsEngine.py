#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

import sys, os, glob, optparse, re, shutil, subprocess, string, time

import json, urllib2, requests, argparse

import apiclient
from apiclient import discovery
import httplib2
from oauth2client import client
from oauth2client import file as oauth2client_file
from oauth2client import tools

# Authorization codes

#TABLE_ID      = 'YOUR_TABLE_ID'



def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Tool for uploading raster images to Google Maps Engine
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg



#--------------------------------------------------------------------------------

def getCredentials(redo=False):

    # Parse OAuth command line arguments
    parser = argparse.ArgumentParser(parents=[tools.argparser])
    flags = parser.parse_args()
    
    print 'Looking for cached credentials'
    # Check if we already have a file with OAuth credentials
    storage = oauth2client_file.Storage('mapsengine.dat')
    credentials = storage.get()
    if credentials is None or credentials.invalid or redo:
        print 'Getting credentials'
        # Start local server, redirect user to authentication page, receive OAuth
        # credentials on the local server, and store credentials in file
        flow = client.OAuth2WebServerFlow(
            client_id=CLIENT_ID,
            client_secret=CLIENT_SECRET,
            scope='https://www.googleapis.com/auth/mapsengine',
            user_agent='Google-MapsEngineApiSample/1.0')
        credentials = tools.run_flow(flow, storage, flags)
      
    return credentials

def authorize(redo=False):

    credentials = getCredentials(redo)    
    
    # Set up discovery with authorized credentials
    http = httplib2.Http()
    http = credentials.authorize(http)
    
    service = apiclient.discovery.build('mapsengine', 'v1', http=http, developerKey=API_KEY)

    #print dir(service.maps())

    #try:
    #    f = service.maps().list()
    #    print '====> ' + str(f)
    #except apiclient.errors.HttpError, error:
    #  if apiclient.error.resp.status == 401:
    #    print 'Credentials have been revoked, re-obtaining them.'
    
    # TODO: Come up with a check that works so we can force new credentials!    
        #credentials = getCredentials(True)    
        #http        = httplib2.Http()
        #http        = credentials.authorize(http)       
        #service     = apiclient.discovery.build('mapsengine', 'v1', http=http, developerKey=API_KEY)

    jsonCredentials = credentials.to_json()
    jsonDict = json.loads(jsonCredentials)
    #for key, value in jsonDict.iteritems():
    #    print key + ' ---> ' + str(value)
    
    bearerToken = jsonDict['access_token']
    #print bearerToken
    return bearerToken
    

## Read the location of every Feature in draft version of Table.
#features = service.tables().features()
#request = features.list(id=TABLE_ID)
#while request is not None:
#  resource = request.execute()
#  for feature in resource['features']:
#    print feature['geometry']['coordinates']
#
#  # Is there an additional page of features to load?
#  request = features.list_next(request, resource)

def createRasterAsset(bearerToken, inputFile): #TODO: Pass in a file list?
    
    url = 'https://www.googleapis.com/mapsengine/v1/rasters/upload'
    
    justFilename = os.path.basename(inputFile)
    data = ( 
    {
      #"projectId": "298099604529",  # REQUIRED, taken from Google API dashboard
      "projectId": "04070367133797133737",  # REQUIRED, taken from Maps Engine URL
      "name": "TEST",  # REQUIRED
      "description": "TODO",
      "files": [ # REQUIRED
        { "filename": justFilename }
      ],
      #"acquisitionTime": {
      #  "start": "2010-01-01T12:00:00Z",
      #  "end": "2010-12-01T12:00:00Z",
      #  "precision": "second"
      #},
      "draftAccessList": "Map Editors", # REQUIRED
      "attribution": "NASA Public Domain" # REQUIRED
      #"tags": ["a", "b"],
      #"maskType": "autoMask"
    } )
    print(data)
    tokenString = 'Bearer '+bearerToken
    headers = {'Authorization': tokenString,
               'Content-Type': 'application/json'}
    print(headers)
    
    response = requests.post(url, data=json.dumps(data), headers=headers)
    
    print response.text
    
    print('Received status code ' + str(response.status_code))
    if response.status_code == 401:
        print('Error: Unauthorized access!')
    if response.status_code == 403:
        print('Error: Forbidden operation!')
    if response.status_code != 200:
        return (False, response.status_code)
    
    jsonDict = json.loads(response.text)
    #for key, value in jsonDict.iteritems():
    #    print key + ' ---> ' + str(value)
       
    # Return the asset ID from the response
    return (True, jsonDict['id'])
    


def uploadFile(bearerToken, assetId, filename):

    # Check input image
    if not os.path.exists(filename):
        raise Exception('Input image file ' + filename + ' is missing!')
    imageSizeBytes = os.path.getsize(filename)

    # Set up POST request
    justFilename = os.path.basename(filename)
    url = 'https://www.googleapis.com/upload/mapsengine/v1/rasters/'+str(assetId)+'/files?filename='+justFilename
    tokenString = 'Bearer '+bearerToken
    headers = {'Authorization': tokenString,
               'Content-Type': 'image/tiff',
               'Content-Length': str(imageSizeBytes)}

    print headers
    print url
    print filename

    fileList = {'file': open(filename, 'rb')}

    response = requests.post(url, headers=headers, files=fileList)
    
    print response.text
    
    # Check response status code
    print('Received status code ' + str(response.status_code))
    if response.status_code != 204:
        print 'Failed to upload file!'
        return False
    print 'File upload started successfully!'
    return True
 

    # To check upload progress:
    # GET https://www.googleapis.com/mapsengine/v1/rasters/{raster_ID}
    # Authorization: Bearer {token}


def main():

    print ('#################################################################################')
    print ("Running GoogleMapsEngine.py")

    #try:
    #try:
    usage = "usage: GoogleMapsEngine.py <input image> [--keep][--manual]\n  "
    parser = optparse.OptionParser(usage=usage)



    parser.add_option("--manual", action="callback", callback=man,
                      help="Read the manual.")
    (options, args) = parser.parse_args()
    
    #if len(args) < 1:
    #    parser.error('Missing required input!')
    #options.inputPath  = args[0]
    #
    #except(optparse.OptionError, msg):
    #    raise Usage(msg)
    
    #options.inputPath = 'means.png'
    options.inputPath = '/home/smcmich1/data/production/NAC_DTM_M151318807_M181974094/results/output-DEM.tif'
    if not os.path.exists(options.inputPath):
        raise Exception('Input file does not exist!')

    startTime = time.time()

    # Get server authorization
    bearerToken = authorize()
    
    # Create empty raster asset request
    (success, assetId) = createRasterAsset(bearerToken, options.inputPath)
    if not success:
        print 'Refreshing access token...'
        bearerToken = authorize(True)
        (success, assetId) = createRasterAsset(bearerToken, options.inputPath)
    if not success:
        raise Exception('Could not get access token!')
    
    if success:
        print 'Created asset ID ' + str(assetId)
    
        # Load a file associated with the asset
        uploadFile(bearerToken, assetId, options.inputPath)





    endTime = time.time()

    print("Finished in " + str(endTime - startTime) + " seconds.")
    print('#################################################################################')
    return 0

    #except(Usage, err):
    #    print(err)
    #    #print(>>sys.stderr, err.msg)
    #    return 2

if __name__ == "__main__":
    sys.exit(main())
