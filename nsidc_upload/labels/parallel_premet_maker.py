
import os, sys, argparse, multiprocessing, icebridge_common
import extract_oib_flight_metadata

'''
Run the premet tool on all the files in a directory in parallel.
- This creates the .premet and .spo files.
'''

SOFTWARE_FOLDER = '/u/smcmich1/icebridge/label_upload_software'

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


def runPremet(camp, platform, tailnum, ortho_folder, name, label_folder):
    '''Create the .premet and .spo file for one input file (ortho and label)'''

    NUM_RETRIES = 3

    premetPath = os.path.join(ortho_folder, name+'.premet')
    spoPath    = os.path.join(ortho_folder, name+'.spo')

    label_name = name.replace('_ortho.tif', '.tif').replace('RDSISCO4',  'RDSISC4')
    premet_label_path = os.path.join(label_folder, label_name+'.premet')
    spo_label_path    = os.path.join(label_folder, label_name+'.spo')

    if os.path.exists(premet_label_path) and os.path.exists(spo_label_path):
        return # Already done with everything

    # The command looks like this:
    full_ortho_path = os.path.join(ortho_folder, name)
    command = ('python %s/RDSISC4_Parser.py "%s" "%s" "%s" %s' %
                   (SOFTWARE_FOLDER, full_ortho_path, camp, platform, tailnum))

    # Try running the command a few times if we did not create all of the files.
    attempt = 0
    while attempt < NUM_RETRIES:

        #print(command)
        os.system(command)
        if (os.path.exists(premetPath) and os.path.exists(spoPath)):
            break # Quit on success
        print 'Failed to generate ' + premetPath + ' and ' + spoPath +', retrying!'
        attempt += 1

    if attempt >= NUM_RETRIES:
        print 'FAILED to generate ' + premetPath + ' and ' + spoPath
        return
    
    #raise Exception('DEBUG')
    # Now create symlinks for these output files in the label folder
    os.system('ln -s ' + premetPath + ' ' + premet_label_path)
    os.system('ln -s ' + spoPath    + ' ' + spo_label_path   )



def do_work(ortho_folder, label_folder, numProcesses):
    '''Do all of the work except for arg parsing'''

    if not os.path.exists(ortho_folder):
        raise Exception('Input folder does not exist!')

    #if ('cameras' in options.folder) or ('summaries' in options.folder):
    #    print('Nothing to be done in the camera or summaries folders!')
    #    return 0

    (campaign, platform, tailnum) = get_flight_info(ortho_folder)

    files = os.listdir(ortho_folder)

    if numProcesses > 1:
        print('Starting processing pool with ' + str(numProcesses) +' processes.')
        pool = multiprocessing.Pool(numProcesses)
        taskHandles = []

    # Run the tool for each .tif file in the folder
    for f in files:

        if os.path.splitext(f)[1] != '.tif':
            continue

        # Add the command to the task pool
        if numProcesses > 1:
            taskHandles.append(pool.apply_async(runPremet, (campaign, platform, tailnum, ortho_folder, f,
                                                            label_folder)))
        else:
            runPremet(campaign, platform, tailnum, ortho_folder, f, label_folder)


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

        options = parser.parse_args(argsIn)

    except argparse.ArgumentError, msg:
        parser.error(msg)

    do_work(options.ortho_folder, options.label_folder, options.numProcesses)


# Run main function if called from shell
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
