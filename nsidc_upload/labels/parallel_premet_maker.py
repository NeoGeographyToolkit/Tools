
import os, sys, argparse, multiprocessing, icebridge_common
import extract_oib_flight_metadata

'''
Run the premet tool on all the files in a directory in parallel.
- This creates the .premet and .spo files.
'''

SOFTWARE_FOLDER = '/u/smcmich1/icebridge/label_upload_software'

def fileIsNonZero(path):  
    '''Return true if the file exists and is non-empty'''
    if os.path.isfile(path) and (os.path.getsize(path) > 0):
        return True
    else:
        return False

def get_flight_info(folder):

    try:
        name  = None
        parts = folder.split('/')
        for p in parts:
            if 'GR_' in p:
                name = p
                break
        if not name:
            raise Exception('')
    except:
        raise Exception('Failed to parse folder: ' + folder)

    # Look up the flight info
    try:
        (campaign, platform, tailnum) = extract_oib_flight_metadata.main([name])
        return (campaign, platform, tailnum)
    except Exception as e:
        print str(e)
        raise Exception("THE CAMPAIGN, PLATFORM AND TAIL NUMBER COULD NOT BE FOUND FOR " + name)


def runPremet(camp, platform, tailnum, ortho_folder, label_folder, kml_file, name):
    '''Create the .premet and .spo file for one input file (ortho and label)'''

    NUM_RETRIES = 3

    parts = name.split('_')
    if len(parts) < 6:
        print('Skipping invalid file name: ' + name)
        return
    frame = int(parts[4])
    #label_name = name.replace('_ortho.tif', '.tif').replace('RDSISCO4',  'RDSISC4')
    ortho_name = name.replace('.tif', '_ortho.tif').replace('RDSISC4',  'RDSISCO4')
    full_label_path = os.path.join(label_folder, name)
    premet_label_path = full_label_path + '.premet'
    spo_label_path    = full_label_path + '.spo'

    if fileIsNonZero(premet_label_path) and fileIsNonZero(spo_label_path):
        return # Already done with everything

    if ortho_folder is not None:
        full_ortho_path = os.path.join(ortho_folder, ortho_name)

        # If there is no matching ortho name for this file, try to find a higher frame
        #  number ortho file that we can steal the position info from.
        orig_frame = frame
        inc = 1
        while not os.path.exists(full_ortho_path):
            frame += inc
            if frame > 30000: # Try searching downward from the first frame this time
                inc = -1
                frame = orig_frame - 1

            if frame < 0: # Tried both directions, hopefully this never happens!
                print('ERROR: Failed to generate premet for file: ' + full_label_path)
                return

            parts[4] = ('%05d' % (frame))
            new_name = '_'.join(parts)
            ortho_name = new_name.replace('.tif', '_ortho.tif').replace('RDSISC4',  'RDSISCO4')
            print('Trying to copy nav data from file: ' + ortho_name)
            full_ortho_path = os.path.join(ortho_folder, ortho_name)

            #if not os.path.exists(full_ortho_path):  # HACK only for flights that span a day and are missing orthos!
            #    ortho_name = ortho_name.replace('_12_', '_11_') # There are actually a bunch of these flights!
            #    full_ortho_path = os.path.join(ortho_folder, ortho_name)
            print('full_ortho_path = ' + full_ortho_path)

        premet_ortho_path = full_ortho_path + '.premet'
        spo_ortho_path    = full_ortho_path + '.spo'

        # The command looks like this:
        command = ('python %s/RDSISC4_Parser.py "%s" "%s" "%s" %s' %
                       (SOFTWARE_FOLDER, full_ortho_path, camp, platform, tailnum))

    else: # Need to use the KML file here
        command = ('python %s/RDSISC4_Parser.py "%s" "%s" "%s" %s' %
                       (SOFTWARE_FOLDER, full_label_path, camp, platform, tailnum))


    print(command)

    # Try running the command a few times if we did not create all of the files.
    attempt = 0
    while attempt < NUM_RETRIES:

        os.system(command)

        if ortho_folder is None:
            # Use the KML parsing script to generate the .spo file
            cmd2 = ('%s/ParseKmlSingle.sh %s %s %s' % (SOFTWARE_FOLDER, full_label_path, name, kml_file))
            print(cmd2)
            os.system(cmd2)

            if (fileIsNonZero(premet_label_path) and fileIsNonZero(spo_label_path)):
                break # Quit on success
            print 'Failed to generate ' + premet_label_path + ' and ' + spo_label_path +', retrying!'
        else:
            if (fileIsNonZero(premet_ortho_path) and fileIsNonZero(spo_ortho_path)):
                break # Quit on success
            print 'Failed to generate ' + premet_ortho_path + ' and ' + spo_ortho_path +', retrying!'
        attempt += 1

    if attempt >= NUM_RETRIES:
        print 'FAILED to generate output files!'
        return
    
    #raise Exception('DEBUG')
    if ortho_folder is not None:
        # Now create symlinks for these output files in the label folder
        os.system('ln -s ' + premet_ortho_path + ' ' + premet_label_path)
        os.system('ln -s ' + spo_ortho_path    + ' ' + spo_label_path   )



def do_work(ortho_folder, label_folder, kml_file, numProcesses):
    '''Do all of the work except for arg parsing'''

    if (not ortho_folder) and (not kml_file):
        raise Exception('Need either the ortho folder or the kml file!')

    #if ('cameras' in options.folder) or ('summaries' in options.folder):
    #    print('Nothing to be done in the camera or summaries folders!')
    #    return 0

    (campaign, platform, tailnum) = get_flight_info(label_folder)

    if ortho_folder is not None:
        ortho_files = os.listdir(ortho_folder)
    else:
        ortho_files = []
    label_files = os.listdir(label_folder)
    print(str(len(ortho_files)) + ' files in ortho folder.')
    print(str(len(label_files)) + ' files in label folder.')

    if numProcesses > 1:
        print('Starting processing pool with ' + str(numProcesses) +' processes.')
        pool = multiprocessing.Pool(numProcesses)
        taskHandles = []

    # Run the tool for each .tif file in the folder
    for f in label_files:

        if os.path.splitext(f)[1] != '.tif':
            continue

        # Add the command to the task pool
        if numProcesses > 1:
            taskHandles.append(pool.apply_async(runPremet, (campaign, platform, tailnum, ortho_folder,
                                                            label_folder, kml_file, f)))
        else:
            runPremet(campaign, platform, tailnum, ortho_folder, label_folder, kml_file, f)


    if numProcesses > 1:
        # Wait for all the tasks to complete
        print('Finished adding ' + str(len(taskHandles)) + ' tasks to the pool.')
        icebridge_common.waitForTaskCompletionOrKeypress(taskHandles, interactive=False)

        # All tasks should be finished, clean up the processing pool
        icebridge_common.stopTaskPool(pool)

    print('Jobs finished.')


def main(argsIn):

    try:
        usage = '''usage: parallel_premet.py'''

        parser = argparse.ArgumentParser(usage=usage)

        # Data selection optios
        parser.add_argument('--ortho-folder', dest='ortho_folder',
                            help="Folder with ortho images.")

        parser.add_argument('--label-folder', dest='label_folder',
                            help="Folder with label images.")

        parser.add_argument('--num-processes', dest='numProcesses', type=int,
                            default=4,
                            help='How many processes to start at the same time.')

        parser.add_argument('--kml-file', dest='kml_file',
                            help="KML file with frame locations.", default=None)


        options = parser.parse_args(argsIn)

    except argparse.ArgumentError, msg:
        parser.error(msg)

    do_work(options.ortho_folder, options.label_folder, options.kml_file, options.numProcesses)


# Run main function if called from shell
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
