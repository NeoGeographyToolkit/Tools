
import os, sys, logging, datetime, multiprocessing, icebridge_common, shutil

import threading, time

import flight_list # Contains AN_FLIGHTS and GR_FLIGHTS

#===============================================================================
# Constants

COMPLETED_DATES_FILE = '/u/smcmich1/icebridge/upload_software/completed_dates.txt'
FAILED_DATES_FILE    = '/u/smcmich1/icebridge/upload_software/failed_dates.txt'
SOFTWARE_FOLDER      = '/u/smcmich1/icebridge/upload_software/'
SIPS_LOG_FOLDER      = os.path.join(SOFTWARE_FOLDER, 'logs')

INITIAL_UNPACK_FOLDER = '/nobackup/smcmich1/icebridge/nsidc_uploads/temp_tar_unpack'

ASP_PATH = '/u/smcmich1/programs/StereoPipeline-2.6.0-2018-02-06-x86_64-Linux/bin/'

EMAIL_ADDRESS = 'scott.t.mcmichael@nasa.gov'

# Controls access to the two DATES files.
FLIGHT_DATES_RW_LOCK = threading.Lock()

#===============================================================================
# Classes

class Date:

    def __init__(self, year, month, day):
        self.year  = int(year)
        self.month = int(month)
        self.day   = int(day)

    def __str__(self):
        return ('%04d%02d%02d' % (self.year, self.month, self.day))

    def yyyymmdd(self):
        return ('%04d%02d%02d' % (self.year, self.month, self.day))

    def yyyy_mm_dd(self):
        return ('%04d_%02d_%02d' % (self.year, self.month, self.day))

    def yyyy_mm(self):
        return ('%04d_%02d' % (self.year, self.month))

class Campaign:

    def __init__(self, site, year):
        self.site = site
        self.year = year

    def getAllDates(self):
        '''Return a list of all the valid dates for this campaign'''

        if self.site == 'AN':
            input_list = flight_list.AN_FLIGHTS
        else:
            input_list = flight_list.GR_FLIGHTS
        flights = [f for f in input_list if f.startswith(self.year)]

        dates = []
        for f in flights:
            year  = f[0:4]
            month = f[4:6]
            day   = f[6:8]
            dates.append(Date(year,month,day))

        return dates

    def __str__(self):
        '''Return the string representation'''
        return self.site + '_' + self.year



def upload_log_and_cleanup(dem_folder, ortho_folder, dem_summary_folder, ortho_summary_folder,
                           unpack_prefix, has_dem_summary, has_ortho_summary, 
                           dem_tarball, ortho_tarball, camera_tarball, summary_tarball,
                           remote_folder, date, logger):
    '''Called by the worker thread in UploadManager'''

    print 'Ready to upload folder ' + dem_folder
    print 'Ready to upload folder ' + ortho_folder
    #return # DEBUG!!!!

    try:
        # Upload the data
        logger.info('Beginning data upload for flight ' + str(date))
        uploadFolderToNSIDC(dem_folder,   'dem/'  +remote_folder, logger)
        uploadFolderToNSIDC(ortho_folder, 'ortho/'+remote_folder, logger)
        if has_dem_summary:
            uploadFolderToNSIDC(dem_summary_folder,   'dem/'  +remote_folder, logger)
        if has_ortho_summary:
            uploadFolderToNSIDC(ortho_summary_folder, 'ortho/'+remote_folder, logger)
        logger.info('Data upload finished for flight ' + str(date))

        success = has_dem_summary and has_ortho_summary

    except Exception as e:
        success = False
        logger.error('Caught exception for date ' + str(date) +'\n' + str(e))

    if success:

        # Log the input tarballs we used and whether we had all summary files.
        updateLogFile(COMPLETED_DATES_FILE, date, dem_tarball, ortho_tarball, 
                      camera_tarball, summary_tarball,
                      has_dem_summary, has_ortho_summary)
        subject = 'COMPLETED flight date: ' + str(date)
        logger.info(subject)

    else:
        updateLogFile(FAILED_DATES_FILE, date, dem_tarball, ortho_tarball, 
                      camera_tarball, summary_tarball,
                      has_dem_summary, has_ortho_summary)
        subject = 'FAILED to process flight date '+ str(date)
        logger.error(subject)

    sendEmail(EMAIL_ADDRESS, subject, 'NT')

    # Clean up the temporary folders
    #raise Exception('DEBUG')

    if success:
        logger.info('Ready to delete folders: ' + unpack_prefix)
        cmd = 'rm -rf ' + unpack_prefix + '*'
        logger.info(cmd)
        os.system(cmd)




class UploadManager():
    '''Class to keep uploading data in the background while the main process starts a new flight.'''


    def __init__(self):
        self._worker = None

    def __del__(self):
        self.cleanup()

    def cleanup(self):
        if self._worker != None:
            self._worker.join()
        self.worker = None

    def uploadFlight(self, dem_folder, ortho_folder, dem_summary_folder, ortho_summary_folder,
                     unpack_prefix, has_dem_summary, has_ortho_summary, 
                     dem_tarball, ortho_tarball, camera_tarball, summary_tarball,
                     remote_folder, date, logger):
        '''Upload the flight in a separate thread.  If another flight is still being uploaded,
           blocks until that upload is finished.'''

        # Block here until we are not busy with another upload.
        if self._worker != None:
            self._worker.join()

        # Set up a working thread with the information
        self._worker = threading.Thread(target=upload_log_and_cleanup,
                                        args=(dem_folder, ortho_folder, dem_summary_folder,
                                              ortho_summary_folder, unpack_prefix,
                                              has_dem_summary, has_ortho_summary,
                                              dem_tarball, ortho_tarball, 
                                              camera_tarball, summary_tarball,
                                              remote_folder, date, logger))
        # Let the worker thread run on its own
        self._worker.start()
               
        return




#===============================================================================
# Functions

def sendEmail(address, subject, body):
    '''Send a simple email from the command line'''
    # Remove any quotes, as that confuses the command line.
    subject = subject.replace("\"", "")
    body    = body.replace("\"", "")
    try:
        cmd = 'mail -s "' + subject + '" ' + address + ' <<< "' + body + '"'
        os.system(cmd)
    except Exception as e:
        print("Could not send mail.")


def getLatestTarFileMatchingDate(dirs, date):
    '''Find the most recent tar file containing a date in the given folders'''
    date_string = str(date)
    candidates = []
    for d in dirs:
        # Get all matching tar files (possibly multiple versions)
        tars = os.listdir(d)
        new_candidates = [f for f in tars if ('.tar' in f) and (date_string in f)]
        candidates = candidates + [os.path.join(d, f) for f in new_candidates]

        # Ignore files manually marked not to use!
        candidates = [c for c in candidates if 'old' not in c]

    if not candidates:
        raise Exception('No tarballs found for date ' + str(date))

    # The last file alphabetically is the latest one
    return sorted(candidates)[-1]

# Functions to find needed tarballs.
def findDemTarball(campaign, date):
    dirs = ['/u/smcmich1/icebridge/output', '/u/smcmich1/icebridge/oleg_dems']
    return getLatestTarFileMatchingDate(dirs, date)

def findOrthoTarball(campaign, date):
    dirs = ['/u/smcmich1/icebridge/ortho', '/u/smcmich1/icebridge/oleg_ortho']
    return getLatestTarFileMatchingDate(dirs, date)

def findCameraTarball(campaign, date):
    dirs = ['/u/smcmich1/icebridge/camera', '/u/smcmich1/icebridge/oleg_cameras']
    return getLatestTarFileMatchingDate(dirs, date)

def findSummaryTarball(campaign, date):
    dirs = ['/u/smcmich1/icebridge/summaries', '/u/smcmich1/icebridge/oleg_summaries']
    return getLatestTarFileMatchingDate(dirs, date)

def fetchTarballsFromTapes(campaign, date_list, logger):
    '''Request that all of the tarballs we will need for this run be loaded from tape.'''

    logger.info('Locating all the tarballs needed for ' + str(len(date_list)) + ' dates.')

    # Find all the tarballs we will need
    needed_files = []
    for date in date_list:
        try:
            dem_tarball     = findDemTarball    (campaign, date)
            ortho_tarball   = findOrthoTarball  (campaign, date)
            camera_tarball  = findCameraTarball (campaign, date)
            summary_tarball = findSummaryTarball(campaign, date)
            needed_files.append(dem_tarball)
            needed_files.append(ortho_tarball)
            needed_files.append(camera_tarball)
            needed_files.append(summary_tarball)
        except:
            logger.error('Error finding all tarballs for date: ' + str(date))

    logger.info('Requesting that these dates be loaded from tape!')

    # Build a command to fetch them all at once.
    cmd = 'dmget '
    for f in needed_files:
        cmd += f + ' '
    cmd += '&' # Run this command in the background so we can start processing as soon as files are ready.
    logger.info(cmd)
    os.system(cmd)


def unpackTarAndGetFileList(tarPath, storage_folder, flight_title, logger, isSummary=False):
    '''Extract the tif files from a tarball into a specified folder.'''

    logger.info('Unpacking tar file: ' + tarPath)

    if os.path.exists(storage_folder):
        logger.info('Storage folder already exists, skipping unpack.')
    else:
        # Each flight uses a different temp unpack location
        this_unpack_folder = os.path.join(INITIAL_UNPACK_FOLDER, flight_title)
        os.system('mkdir -p ' + this_unpack_folder)
        cmd = 'tar -xf ' + tarPath + ' --directory ' + this_unpack_folder
        print cmd
        logger.info(cmd)
#        os.system(cmd)

        logger.info('Finished tar unpack command, looking for output...')

        possible_directories = ['tarAssembly', 'processed', 'camera', 'summary', flight_title]
        file_dir   = []
        top_folder = os.path.join(this_unpack_folder, flight_title)
        for d in possible_directories:
            test_dir = os.path.join(top_folder, d)
            print(test_dir)
            if os.path.exists(test_dir):
                file_dir = test_dir
                break
            test_dir = os.path.join(this_unpack_folder, d)
            print(test_dir)
            if os.path.exists(test_dir):
                file_dir = test_dir
                break

        if not file_dir:
            raise Exception('ERROR: Did not find unpack folders for storage folder ' + storage_folder)

        logger.info('Found data in: ' + file_dir + ', moving to ' + storage_folder)

        # Move all the data files into a new directory
        cmd = 'mv ' + file_dir +' '+ storage_folder
        print cmd
        logger.info(cmd)
        os.system(cmd)

        # Delete the unpack folder.
        cmd = 'rm -rf ' + this_unpack_folder
        print cmd
        logger.info(cmd)
        os.system(cmd)

    logger.info('Retrieving the file list...')

    # Get all the .tif files in the folder
    # - Also need to handle cases where files are in a subfolder.
    #   In those cases we record the subfolder path and will get the file later.
    all_file_list = os.listdir(storage_folder)
    file_list     = []
    bad_file_list = []
    needed_file_types = ['.tif', '.tsai', '.jpg', '.jpeg']
    for f in all_file_list:
        full_path = os.path.join(storage_folder, f)
        ext       = os.path.splitext(f)[1]
        if (ext in needed_file_types) or (os.path.isdir(full_path) and len(f) > 3):
            # Empty image files are a problem, but empty camera files
            #  are ok since we only use them for the file names.
            if (os.path.getsize(full_path) == 0) and (ext != '.tsai'):
                if isSummary: # We will just regenerate these later
                    print 'Deleting empty summary file: ' + f
                    os.remove(full_path)
                else:
                    bad_file_list.append(full_path)
                    print('After unpack, got empty file: ' + f) 
            else: # A good file!
                file_list.append(full_path)

    num_bad_files = len(bad_file_list)
    logger.info('Num bad files = ' + str(num_bad_files))
    if num_bad_files > 0:
        raise Exception('Quitting because of missing files after unpacking ' + tarPath)
    return file_list


def add_timestamps_to_files(input_files, camera_files, postfix, browse):
    '''Update the input file names to include timestamps'''

    if not input_files:
        return

    # Associate each camera file with its frame number
    camera_frames = {}
    for c in camera_files:
        parts = os.path.basename(c)[0:-5].split('_')
        frame = int(parts[3])
        camera_frames[frame] = parts


    # Rename each of the DEM files to this format:
    #   IODEM3_20091016_17534868_02172_DEM.tif
    # Ortho files start with IODIM3.
    prefix = 'IODEM3_'
    if postfix == 'ORTHO':
         prefix = 'IODIM3_'
    input_dir     = os.path.dirname(input_files[0])
    missing_files = False
    for in_file in input_files:

        fname = os.path.basename(in_file)

        if prefix in fname: # Skip already converted files
            continue

        parts  = os.path.splitext(fname)[0].split('_')
        if len(parts) > 1:
            frame  = int(parts[1])
        else: # Handle old files with just the frame number
            frame  = int(parts[0])

        try:
            cam_parts = camera_frames[frame]
        except KeyError:
            print('Missing camera file for input image: ' + in_file)
            missing_files = True
            continue
        new_name = (prefix + cam_parts[1] +'_' + cam_parts[2]
                     +'_' + cam_parts[3] +'_'+ postfix)
        if browse:
            if postfix == 'ORTHO':
                new_name += '.jpg'
            else: # DEM
                new_name += '_browse.tif'
        else:
            new_name += '.tif'
        new_path = os.path.join(input_dir, new_name)

        # Sometimes the file is inside a folder where the folder has the frame name
        if os.path.isdir(in_file):
            sub_files = os.listdir(in_file)
            if len(sub_files) != 1:
                raise Exception('Too many subfiles for folder: ' + in_file)
            cmd = 'mv ' + os.path.join(in_file, sub_files[0]) +' '+ new_path
            #print cmd
            os.system(cmd)
            cmd = 'rm -rf ' + in_file # Clean up the empty folder
            #print cmd
            os.system(cmd)

        else: # Normal file
            cmd = 'mv ' + in_file +' '+ new_path
            #print cmd
            os.system(cmd)

    if missing_files:
        raise Exception('Missing at least one camera file, check the log!')

def makeConfigCopy(source_path, output_path, data_folder, browseFolder, isDim):
    '''Make a copy of the the .cfg file with the data folder inserted.'''

    os.system('rm -f ' + output_path) # Always start fresh!

    source_file = open(source_path, 'r')
    output_file = open(output_path, 'a')

    # Point to the correct config files
    mcfName = 'IODEM3#001.MCF'
    pushName = 'push_DEM.cfg'
    if isDim:
        mcfName  = 'IODIM3#001.MCF'
        pushName = 'push_DIM.cfg'

    for line in source_file:
        line = line.replace('MCF_REPLACE',  mcfName)
        line = line.replace('PUSH_REPLACE', pushName)
        if browseFolder:
            line = line.replace('BROWSE_PRESENT_REPLACE', 'Y')
        else:
            line = line.replace('BROWSE_PRESENT_REPLACE', 'N')
        line = line.replace('BROWSE_FOLDER_REPLACE', browseFolder + '/')
        if isDim: # ORTHO
            line = line.replace('BROWSE_EXT_REPLACE', '.jpg')
        else: # DEM
            if browseFolder:
                line = line.replace('BROWSE_EXT_REPLACE', '_browse.tif')
            else: # Original file
                line = line.replace('BROWSE_EXT_REPLACE', '.tif')
        line = line.replace('REPLACE',     data_folder)
        output_file.write(line)
    source_file.close()
    output_file.close()

def createSummaryFile(inFile, outFile, isOrtho):
   '''Create one summary file'''

   if isOrtho:
       cmd = ASP_PATH + 'gdal_translate -scale -outsize 25% 25% -of jpeg ' + inFile +' '+ outFile
       print cmd
       os.system(cmd)
   else:
       # Hillshade then downsample.
       tempPath = outFile + '_temp.tif'
       cmd = ASP_PATH + 'hillshade ' + inFile +' -o ' + tempPath
       print cmd
       os.system(cmd)
       cmd = (ASP_PATH +'gdal_translate '+ tempPath +' '+ outFile+
              ' -of GTiff -outsize 40% 40% -b 1 -co "COMPRESS=JPEG"')
       print cmd
       os.system(cmd)
       os.system('rm ' + tempPath)

def replaceMissingSummaryFiles(folder, summaryFolder, isOrtho, logger):
    '''Recreate any missing summary files'''

    inFiles      = os.listdir(folder)
    summaryFiles = os.listdir(summaryFolder)

    # Associate each summary file with its frame number
    summaryFrames = {}
    for c in summaryFiles:
        parts = os.path.basename(c).split('_')
        frame = int(parts[3])
        summaryFrames[frame] = c

    # See if any input files are missing their summary frame
    missingSummaryList = []
    for f in inFiles:
        # Skip metadata files
        ext = os.path.splitext(f)[1]
        if ext != '.tif':
            continue

        parts = os.path.basename(f).split('_')
        frame = int(parts[3])

        success = False
        try:
            summary  = summaryFrames[frame]
            # File exists, make sure it is not empty.
            fullPath = os.path.join(summaryFolder, summary)
            if (os.path.getsize(fullPath) > 0):
                success = True
                continue
        except KeyError: # Other errors are unexpected
            pass # File does not exist

        if not success:
            logger.info('Detected missing or empty summary file for: ' + f)

            # Get the desired summary file name
            if isOrtho:
                summaryName = f.replace('.tif', '.jpg')
            else: # DEM
                summaryName = f.replace('.tif', '_browse.tif')
            summaryPath = os.path.join(summaryFolder, summaryName)
            inputPath   = os.path.join(folder, f)
            missingSummaryList.append((inputPath, summaryPath))

    # If we are missing too many files fail out, otherwise regenerate them.

    MAX_RECREATE_FILES = 200
    numMissingFiles = len(missingSummaryList)
    if numMissingFiles > MAX_RECREATE_FILES:
        raise Exception('Missing too many summary files: ' + str(numMissingFiles))
    else:
        logger.info('Will regenerate ' + str(numMissingFiles) + ' missing summary files.')

    for (inputPath, summaryPath) in missingSummaryList:

        # Recreate the missing file.
        createSummaryFile(inputPath, summaryPath, isOrtho)
        if os.path.exists(summaryPath) and (os.path.getsize(summaryPath) > 0):
            logger.info('Successfully created summary file: ' + summaryPath)
        else:
            logger.error('Failed to create summary file: ' + summaryPath)
            raise Exception('Failed to create summary file: ' + summaryPath)


def verifyMetadataFiles(folder, extensions, logger):
   '''Raise an exception if all of the required metadata files for upload are not present'''

   new_files = os.listdir(folder)
   counts = {}
   counts['tif'] = 0
   for e in extensions:
       counts[e] = 0

   for f in new_files: # Count up all the file extensions
       (base, ext) = os.path.splitext(f)
       base = os.path.join(folder, base)

       if ext == '.tif':
           counts['tif'] += 1

           # For each TIF file, log if the other files are missing so
           #  we can find the bad ones.
           for e in extensions:
               ePath = base + '.tif.' + e

               if not os.path.exists(ePath):
                   logger.error('Missing output file: ' + ePath)

       # Accumulate the proper extension
       for e in extensions:
           if ext == '.'+e:
               counts[e] += 1

   # Make sure all the counts are the same
   for e in extensions:
       if counts['tif'] != counts[e]:
           msg = ('Error: in folder ' + folder + ', counts = ' + str(counts))
           logger.error(msg)
           raise RuntimeError(msg)


def runMetGen(configPath):
    '''Used to run this in a Process'''
    cmd = SOFTWARE_FOLDER + 'MetGen.Start -config ' + configPath
    print cmd
    os.system(cmd)


def uploadFolderToNSIDC(folder, remote_folder, logger):
    '''Send the data to NSIDC!'''

    # Push the directory to NSIDC
    remoteDirPath = os.path.join('/incoming', 'Ames', remote_folder)
    logger.info("Storing at NSIDC in: " + remoteDirPath)

    # Not using remote delete since we are pushing two local folders
    #  into a single remote folder.

    #cmd = ('lftp -e "mirror -P 20 -c -R -vvv --delete --delete-first ' + folder +
    #       ' ' + remoteDirPath + ' -i \'\.(tif|jpg|PDR|met)$\'; bye\" -u ' + auth)
    cmd = ('lftp -e "mirror -P 20 -c -R -vvv ' + folder +
           ' ' + remoteDirPath + ' -i \'\.(tif|jpg|PDR|met)$\'; bye\" -u ' + auth)
    logger.info(cmd)

    #raise Exception('DEBUG')

    status = os.system(cmd)
    logger.info("LFTP status: " + str(status))
    if status != 0:
        raise Exception("Problem pushing folder: " + folder)

def updateLogFile(log_file, date, dem_tarball, ortho_tarball, camera_tarball, summary_tarball,
                  dem_summary_present, ortho_summary_present):
    '''Add an entry to the log file if it is not already there'''
    if checkLogFile(log_file, date):
        return

    with FLIGHT_DATES_RW_LOCK: # Grab the lock before writing
        with open(log_file, 'a') as f:
            s = ('%s, %s, %s, %s, %s, %s, %s' % 
                  (str(date), dem_tarball, ortho_tarball, camera_tarball, summary_tarball,
                  dem_summary_present, ortho_summary_present))
            f.write(s + '\n')

def checkLogFile(log_file, date):
    '''Return true if the date is in the file'''
    with FLIGHT_DATES_RW_LOCK: # Grab the lock before reading
        with open(log_file, 'r') as f:
            for line in f:
                parts = line.split()
                if not parts:
                    continue
                if str(date) in parts[0]:
                    return True
        return False


#========================================================
def main():

    # Take the location and year as an input.
    if len(sys.argv) != 4:
        print 'Usage: oib_uploader.py <site> <year> <unpack_folder>'
        return -1

    this_campaign = Campaign(sys.argv[1], sys.argv[2])
    unpack_folder = sys.argv[3]

    # Set up logging
    logFolder = '/u/smcmich1/icebridge/logs'
    logLevel = logging.INFO
    logger   = icebridge_common.setUpLogger(logFolder, logLevel, "push_"+str(this_campaign))
    logger.info("Logging in: " + logFolder)

    # TODO: Maybe we want to keep this around in the future.
    # SIPS internally records which files it has already processed and will not reprocess
    #  them even if the output files are missing, so delete it to avoid trouble.
    logger.info('Deleting SIPS log folder.')
    os.system('rm -rf ' + SIPS_LOG_FOLDER)

    valid_dates = this_campaign.getAllDates()

    # Filter out the dates we have already uploaded.
    dates_to_upload = []
    for d in valid_dates:
        if ( (not checkLogFile(COMPLETED_DATES_FILE, d)) and 
             (not checkLogFile(FAILED_DATES_FILE, d))    ):
            dates_to_upload.append(d)
    for d in dates_to_upload:
        print(d)

    if not dates_to_upload:
        print 'No dates to upload!'

    #dates_to_upload = [Date(2010,05,20)] # DEBUG

    # Start retrieving all of the files we will need from tape in the background.
    # - This takes a long time but it is more efficient to do all at once and
    #   we can start processing as soon as required files are finished.
    fetchTarballsFromTapes(this_campaign, dates_to_upload, logger)

    UPLOAD_LIMIT  = 4 # Upload this many flights before stopping.
    num_uploaded  = 0

    # This object holds on to the uploading (as opposed to prep) thread.
    uManager = UploadManager()

    # Loop through all the dates that have not been uploaded yet.
    for date in dates_to_upload:

        logger.info('Uploading date: ' + str(date))

        # Make a temporary folder
        print 'Using unpack folder: ' + unpack_folder
        os.system('mkdir -p ' + unpack_folder)
        os.system('cd ' + unpack_folder)

        dem_summary_files   = [] # Init for later
        ortho_summary_files = []
        dem_tarball     = None
        ortho_tarball   = None
        camera_tarball  = None
        summary_tarball = None

        success = False
        try: # Make one attempt on a flight, we will have to manually check failures.

            # Unpack the DEM and ORTHO tarballs into new folders.
            logger.info('Unpacking tarballs...')
            flight_title = this_campaign.site +'_'+ date.yyyymmdd()
            folder_prefix  = os.path.join(unpack_folder, flight_title)
            dem_folder     = folder_prefix + '_dems'
            ortho_folder   = folder_prefix + '_orthos'
            camera_folder  = folder_prefix + '_cameras'
            summary_folder = folder_prefix + '_summaries' # Temporary folder until the summary files are split
            dem_summary_folder   = folder_prefix + '_dem_summaries'
            ortho_summary_folder = folder_prefix + '_ortho_summaries'


            dem_tarball     = findDemTarball    (this_campaign, date)
            ortho_tarball   = findOrthoTarball  (this_campaign, date)
            camera_tarball  = findCameraTarball (this_campaign, date)
            summary_tarball = findSummaryTarball(this_campaign, date)

            logger.info('Using output folder prefix: ' + folder_prefix)
            logger.info('Found DEM tarball:     ' + dem_tarball   )
            logger.info('Found ORTHO tarball:   ' + ortho_tarball )
            logger.info('Found CAMERA tarball:  ' + camera_tarball)
            logger.info('Found SUMMARY tarball: ' + summary_tarball )

            dem_files     = unpackTarAndGetFileList(dem_tarball,     dem_folder,     flight_title, logger)
            ortho_files   = unpackTarAndGetFileList(ortho_tarball,   ortho_folder,   flight_title, logger)
            summary_files = unpackTarAndGetFileList(summary_tarball, summary_folder, flight_title, logger, isSummary=True)
            logger.info('Found ' + str(len(dem_files)) +' dem files and ' + str(len(ortho_files)) +' ortho files.')

            # Split out the summary files into dem and ortho files and refresh the file list
            logger.info('Splitting out summary files...')
            os.system('mkdir '+ dem_summary_folder  )
            os.system('mkdir '+ ortho_summary_folder)
            for f in os.listdir(summary_folder):
                absPath = os.path.join(summary_folder, f)
                if 'dem' in f:
                    new_path = os.path.join(dem_summary_folder, f)
                    shutil.move(absPath, new_path)
                if 'ortho' in f:
                    new_path = os.path.join(ortho_summary_folder, f)
                    shutil.move(absPath, new_path)
            dem_summary_files   = [os.path.join(dem_summary_folder,   f) for f in os.listdir(dem_summary_folder  )]
            ortho_summary_files = [os.path.join(ortho_summary_folder, f) for f in os.listdir(ortho_summary_folder)]


            # Handle old cases where the archived files don't have timestamps
            # - Does this catch all the cases?
            last_dem_name   = os.path.basename(dem_files[-1]    )
            last_ortho_name = os.path.basename(ortho_files[-1]  )

            logger.info('Updating main file names with timestamps...')
            # Fetch the camera files, then steal the timestamps from them.
            camera_files = unpackTarAndGetFileList(camera_tarball, camera_folder, flight_title, logger)
            logger.info('Found ' + str(len(camera_files)) +' camera files.')

            if (last_dem_name.startswith('F') or last_dem_name.startswith('batch')):
                add_timestamps_to_files(dem_files,  camera_files, 'DEM', False)
            if (last_ortho_name.startswith('F') or last_ortho_name.startswith('batch')):
                add_timestamps_to_files(ortho_files, camera_files, 'ORTHO', False)


            # The summary files never had timestamps
            logger.info('Updating summary file names with timestamps...')
            add_timestamps_to_files(dem_summary_files,   camera_files, 'DEM',   True)
            add_timestamps_to_files(ortho_summary_files, camera_files, 'ORTHO', True)            

            #raise Exception('DEBUG')
            logger.info('Check for missing DEM summaries...')
            replaceMissingSummaryFiles(dem_folder,   dem_summary_folder,   False, logger)
            logger.info('Check for missing ORTHO summaries...')
            replaceMissingSummaryFiles(ortho_folder, ortho_summary_folder, True,  logger)

            logger.info('Executing uploader prep script...')

            # Execute the BatcherIODEM3PremetSpatialMaker.sh script
            # - This will create some supporting files for each DEM file.
            # - Try a few times in case it fails.
            MAX_NUM_TRIES = 3
            numTries = 0
            while numTries < MAX_NUM_TRIES:
                try:

                    numTries += 1

                    cmd = 'python ' + SOFTWARE_FOLDER + 'Parallel_BatcherIODEM3PremetSpatialMaker.py ' + unpack_folder +' '+ date.yyyymmdd()
                    logger.info(cmd)
                    print(cmd)
                    os.system(cmd)

                    # Check that the preliminary files are all there.
                    verifyMetadataFiles(dem_folder,   ['spo', 'premet'], logger)
                    verifyMetadataFiles(ortho_folder, ['spo', 'premet'], logger)

                except RuntimeError as e:
                    logger.error(str(e))
                    if numTries < MAX_NUM_TRIES:
                        logger.info('Trying Premet gen again...')
                        print 'Trying Premet gen again.....'
                    else:
                        raise Exception('Premet gen failed after ' + str(numTries) + ' attempts!')


            #raise Exception('DEBUG')
            print 'Executing MetGen script on DEMS.....'
            logger.info('Executing MetGen script on DEMS...')

            # Need to update "Primary.cfg" for every input folder.
            PRIMARY_CONFIG_PATH    = SOFTWARE_FOLDER + 'Primary.cfg.src'
            TEMP_CONFIG_PATH_DEM   = SOFTWARE_FOLDER + 'Primary_dem.cfg'
            TEMP_CONFIG_PATH_ORTHO = SOFTWARE_FOLDER + 'Primary_ortho.cfg'

            # Generate more metadata files
            logger.info('Launching MetGen.Start commands...')
            makeConfigCopy(PRIMARY_CONFIG_PATH, TEMP_CONFIG_PATH_DEM,   dem_folder,   dem_summary_folder,   isDim=False)
            makeConfigCopy(PRIMARY_CONFIG_PATH, TEMP_CONFIG_PATH_ORTHO, ortho_folder, ortho_summary_folder, isDim=True )

            # This process sometimes randomly fails on a few frames, so re-run it
            #  in case this happens.
            numTries = 0
            while numTries < MAX_NUM_TRIES:
                try:

                    numTries += 1

                    # This is a java tool so we run both instances in parallel.
                    demMetGenProcess   = multiprocessing.Process(target=runMetGen, args=(TEMP_CONFIG_PATH_DEM,))
                    orthoMetGenProcess = multiprocessing.Process(target=runMetGen, args=(TEMP_CONFIG_PATH_ORTHO,))
                    demMetGenProcess.start()
                    orthoMetGenProcess.start()
                    demMetGenProcess.join()
                    orthoMetGenProcess.join()

                    logger.info('MetGen processes are finished!')      

                    # Check that we created all of the metadata files
                    verifyMetadataFiles(dem_folder,   ['PDR', 'met'], logger)
                    verifyMetadataFiles(ortho_folder, ['PDR', 'met'], logger) # Throws RuntimeError on failure.

                    break # Quit the loop on success.

                except RuntimeError as e:
                    logger.error(str(e))
                    if numTries < MAX_NUM_TRIES:
                        logger.info('Trying MetGen again...')
                        print 'Trying MetGen again.....'
                    else:
                        raise Exception('Metgen failed after ' + str(numTries) + ' attempts!')
 

            logger.info('Verified that all metadata files have been generated.')

            #raise Exception('DEBUG')

            # Start up the uploading thread for this flight,
            #  then start preparing the files for the next flight.
            print 'Calling upload manager!'
            remote_folder = str(date.year) +'/'+ this_campaign.site +'_'+ date.yyyy_mm_dd()
            uManager.uploadFlight(dem_folder, ortho_folder, dem_summary_folder,
                                  ortho_summary_folder, folder_prefix, 
                                  dem_summary_files   != [], ortho_summary_files != [], 
                                  dem_tarball, ortho_tarball, camera_tarball, summary_tarball,
                                  remote_folder, date, logger)

        except Exception as e:
            print 'Caught exception ' + str(e)

            if str(e) != "DEBUG":
                updateLogFile(FAILED_DATES_FILE, date, dem_tarball, ortho_tarball, 
                              camera_tarball, summary_tarball,
                              dem_summary_files   != [],
                              ortho_summary_files != [])

            logger.error('FAILED to process flight date '+ str(date))
            sendEmail(EMAIL_ADDRESS, 'Caught exception for flight: ' + str(date), str(e))


        num_uploaded += 1
        if num_uploaded >= UPLOAD_LIMIT:
            print 'Hit the upload limit!'
            print 'Waiting for the uploader to finish...'
            uManager.cleanup()
            return 0


if __name__=='__main__':
    main()
    sendEmail(EMAIL_ADDRESS, 'Finished running OIB uploader.', 'NT')
