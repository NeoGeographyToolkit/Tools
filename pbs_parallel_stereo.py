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

import sys
import os, glob, re, shutil, subprocess, string, time, errno, optparse
import IrgFileFunctions, IrgIsisFunctions, IrgPbsFunctions, IrgSystemFunctions



def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Calls mapproject in parallel with ISIS cameras.
'''
    sys.exit()


#------------------------------------------------------------------------------

def main(argsIn):

    try:
        usage  = "usage: pbs_parallel_stereo.py (same inputs as parallel_stereo plus new options)"
        parser = IrgSystemFunctions.PassThroughOptionParser(usage=usage) # Use parser that ignores unknown options


        parser.add_option("--num-correlation-nodes",    dest="numCorrelationNodes",   type='int', default=1,
                                              help="Number of nodes to use for the two correlation steps")
        parser.add_option("--num-triangulation-nodes",  dest="numTriangulationNodes", type='int', default=1,
                                              help="Number of nodes to use for the triangulation steps")

        parser.add_option('--node-type',  dest='nodeType', default='wes',
                                         help='Type of processing node to request (wes, san, or ivy)')

        parser.add_option('--group-id',  dest='groupId', default='',
                                         help='GID to charge the hours to [REQUIRED]')


        # This call handles all the specific options for this code.
        (options, args) = parser.parse_args(argsIn)

        # 'args' contains everything for parallel_stereo

        # Check the required positional arguments.
        if not options.groupId:
            parser.error("Must input a group ID to charge to!")

        # Any additional arguments need to be forwarded to the mapproject function
        options.extraArgs = args

    except optparse.OptionError, msg:
        raise Usage(msg)

    #startTime = time.time()
    
    # Currently all the outputs are written to the current directory!
    scriptFolder = os.getcwd()
    pbsPath      = os.path.abspath(os.path.join(scriptFolder, 'mainPbsScript.sh'))
    stdLogPaths  = []
    errLogPaths  = []
    scriptCalls  = []
    for i in range(4):
        stdLogPaths.append(os.path.abspath(os.path.join(scriptFolder, 'stdLog'   +str(i)+'.txt')))
        errLogPaths.append(os.path.abspath(os.path.join(scriptFolder, 'errLog'   +str(i)+'.txt')))
        scriptCalls.append(os.path.abspath(os.path.join(scriptFolder, 'pbsScript'+str(i)+'.sh' )))
    
    # Generate the core PBS string
    cpusPerNode = 12
    if options.nodeType == 'san':
        cpusPerNode = 18
    elif options.nodeType == 'ivy':
        cpusPerNode = 24
        
    # TODO: Allow users to input these times!
    stepHours = ['5:00:00', '40:00:00', "30:00:00", "40:00:00"]
    
    corePbsString = ('qsub -q long -W group_list='+ options.groupId+
                     ' -m eb -S /bin/bash -V -j oe -C '+ scriptFolder)
    
    
    # Generate all of the individual PBS calls
    pbsStrings = []
    
    # Preprocessing stage
    pbsStrings.append('subjob1=$( ' + corePbsString + ' -N pbs_stereo1 -l walltime="'+stepHours[0]+'"'
                                    + ' -e '+ errLogPaths[0] +' -o '+ stdLogPaths[0]
                                    + ' -l select='+str(1)+':ncpus='+str(cpusPerNode)+':model='+options.nodeType
                                    + '    --    '+ scriptCalls[0] +')')
                                    
    # Correlation stage
    pbsStrings.append('subjob2=$( ' + corePbsString + '-N pbs_stereo2 -l walltime="'+stepHours[1]+'"'
                                    + ' -e '+ errLogPaths[1] +' -o '+ stdLogPaths[1]
                                    + ' -l select='+str(options.numCorrelationNodes)+':ncpus='+str(cpusPerNode)+':model='+options.nodeType
                                    + ' -W depend=afterok:$subjob1    --    '+ scriptCalls[1] +')')

    # Filtering stage
    pbsStrings.append('subjob3=$( ' + corePbsString + '-N pbs_stereo3 -l walltime="'+stepHours[2]+'"'
                                    + ' -e '+ errLogPaths[2] +' -o '+ stdLogPaths[2]
                                    + ' -l select='+str(1)+':ncpus='+str(cpusPerNode)+':model='+options.nodeType
                                    + ' -W depend=afterok:$subjob2    --    '+ scriptCalls[2] +')')

    # Triangulation stage
    pbsStrings.append(corePbsString + '-N pbs_stereo4 -l walltime="'+stepHours[3]+'"'
                                    + ' -e '+ errLogPaths[3] +' -o '+ stdLogPaths[3]
                                    + ' -l select='+str(options.numTriangulationNodes)+':ncpus='+str(cpusPerNode)+':model='+options.nodeType
                                    + ' -W depend=afterok:$subjob3    --    '+ scriptCalls[3])

    # Set up the command line for parallel_stereo
    commandList   = ['parallel_stereo',  '--nodes-list', '$PBS_NODEFILE']
    commandList   = commandList + options.extraArgs # Append other options
    commandString = IrgSystemFunctions.argListToString(commandList)

    phases = [' --entry-point 0 --stop-point 1', # Init
              ' --entry-point 1 --stop-point 3', # Correlation
              ' --entry-point 3 --stop-point 4', # Filtering
              ' --entry-point 4 --stop-point 6'] # Triangulation
    
    # Generate a set of four script files
    for i in range(4):
        scriptFile = open(scriptCalls[i], 'w')
        scriptFile.write('#!/bin/bash\n\n')
        thisCommandString = commandString + phases[i]
        scriptFile.write(thisCommandString)
        scriptFile.close()
    
        # Set the script file to be executable
        os.system('chmod +x ' + scriptCalls[i])
    
    # Write the PBS script
    scriptFile = open(pbsPath, 'w')
    scriptFile.write('#!/bin/bash\n\n\n')
    scriptFile.write('# '+ commandString) # Show the implemented command in comments
    for i in range(4):
        scriptFile.write('\n\n\n' + pbsStrings[i])
    scriptFile.close()    
    
    
    ## Clean up temporary files
    #if not options.keep:
    #    IrgFileFunctions.removeFolderIfExists(tempFolder)


    #endTime = time.time()
    #
    #print "Finished in " + str(endTime - startTime) + " seconds."

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
