
import os, sys, argparse, multiprocessing, icebridge_common

'''
Run the premet tool on all the files in a directory in parallel.
- This creates the .premet and .spo files.
'''

SOFTWARE_FOLDER = '/u/smcmich1/icebridge/upload_software'

def runCommand(command):
    '''Run one of the commands from the file'''
    print(command)
    os.system(command) 

def runPremet(command):
    '''Create the .premet and .spo file for one input'''

    NUM_RETRIES = 3

    # The command looks like this:
    #command = ('%s/IODEM3PremetSpatialMaker_single.sh "%s" "%s" "%s" %s %s' %
    #               (SOFTWARE_FOLDER, options.camp, options.platform, options.tailnum,
    #                options.folder, f))

    parts  = argString.split(' ')
    folder = parts[5]
    name   = parts[6]
    premetPath = os.path.join(folder, name+'.premet')
    spoPath    = os.path.join(folder, name+'.spo')

    # Try running the command a few times if we did not create all of the files.
    attempt = 0
    while (attempt < NUM_RETRIES):

        os.system(command)
        if (os.path.exists(premetPath) and os.path.exists(spoPath)):
            return # Quit on success
        print 'Failed to generate ' + premetPath + ' and ' + spoPath +', retrying!'
        attempt += 1

    print 'FAILED to generate ' + premetPath + ' and ' + spoPath

def main(argsIn):

    try:
        usage = '''usage: parallel_premet.py'''
   
        parser = argparse.ArgumentParser(usage=usage)

        # Data selection optios
        parser.add_argument('--folder', dest='folder',
                            help="Folder to process.")

        parser.add_argument('--camp', dest='camp',
                            help="camp.")
        parser.add_argument('--platform', dest='platform',
                            help="platform.")
        parser.add_argument('--tailnum', dest='tailnum',
                            help="tailnum.")


        parser.add_argument('--num-processes', dest='numProcesses', type=int,
                            default=4,
                            help='How many processes to start at the same time.')

        options = parser.parse_args(argsIn)
        
    except argparse.ArgumentError, msg:
        parser.error(msg)

    if not os.path.exists(options.folder):
        raise Exception('Input folder does not exist!')

    if ('cameras' in options.folder) or ('summaries' in options.folder):
        print('Nothing to be done in the camera or summaries folders!')
        return 0

    files = os.listdir(options.folder)

    print('Starting processing pool with ' + str(options.numProcesses) +' processes.')
    pool = multiprocessing.Pool(options.numProcesses)
    taskHandles = []

    # Run the tool for each .tif file in the folder
    for f in files:

        if os.path.splitext(f)[1] != '.tif':
            continue

        # TODO: Don't repeat on files?

        # Add the command to the task pool
        command = ('%s/IODEM3PremetSpatialMaker_single.sh "%s" "%s" "%s" %s %s' % 
                   (SOFTWARE_FOLDER, options.camp, options.platform, options.tailnum,
                    options.folder, f))
        #print command
        taskHandles.append(pool.apply_async(runCommand, (command,)))
        

    # Wait for all the tasks to complete
    print('Finished adding ' + str(len(taskHandles)) + ' tasks to the pool.')
    icebridge_common.waitForTaskCompletionOrKeypress(taskHandles, interactive=False)

    # All tasks should be finished, clean up the processing pool
    icebridge_common.stopTaskPool(pool)
    print('Jobs finished.')



# Run main function if called from shell
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
