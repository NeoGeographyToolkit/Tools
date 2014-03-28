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

"""IrgFileFunctions.py - General functions for handling files"""

import sys, os, re, shutil, subprocess, string, time, errno


def createFolder(path):
    """Creates a folder if it does not already exist"""
    if path == '':
        return
    if not os.path.exists(path):
        os.mkdir(path)

def removeIfExists(path):
    """Removes a file if it exists"""
    try:
        os.remove(path)
    except OSError as e: 
        if e.errno != errno.ENOENT: # Continue if the error is "no such file or directory"
            raise # Re-raise the exception if a different error occured

def removeFolderIfExists(directory):
    """Removes a directory and everything in it"""
    try:
        shutil.rmtree(directory)
    except OSError as e: 
        if e.errno != errno.ENOENT: # Continue if the error is "no such file or directory"
            raise # Re-raise the exception if a different error occured

def getFileLineCount(filePath):
    """Counts up the number of lines in a file"""
    f = open(filePath)
    i = 0
    for line in f:
        i = i + 1
    return i

def checkIfToolExists(toolName):
    """Returns true if the system knows about the utility with this name (it is on the PATH)"""

    # Look for the tool using the 'which' command
    cmd = ['which', toolName]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    translateOut, err = p.communicate()
    

    # Check if that command failed to find the file
    failString = 'no ' + toolName + ' in ('
    if translateOut.find(failString) >= 0:
        raise Exception('Missing requested tool ' + toolName)
    else:
        return True


def getLastGitTag(codePath):
    """Returns the last brief git tag of the repository containing the file""" 

    # Get path to git folder
    codeFolder = os.path.dirname(codePath)
    gitFolder  = os.path.join(codeFolder, '.git/')
    
    # Get the tag using a subprocess call
    cmd = ['git', '--git-dir', gitFolder, 'describe', '--abbrev=0']
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    textOutput, err = p.communicate()
    
    # Check for errors
    if (textOutput.find('fatal:') >=0):
        raise Exception('Error: getLastGitTag failed on code path ' + codePath)

    return textOutput

def tarFileList(fileList, outputPath):
    """Creates a tar file containing a list of files with no absolute paths"""

    # This extra set of commands is needed to strip the absolute path name from each stored file
    cmd = 'tar -jcvf ' + outputPath
    for f in fileList:
        cmd = cmd + ' -C ' + os.path.dirname(f) + ' ' + os.path.basename(f)
    print cmd
    os.system(cmd)


#==================================================
# This class implements a variant of OptionParser which ignores unknown options.

# TODO: Replace existing implementation in stereo_utils.py?
#       --> This version does not make the user treat the arguments differently than normal.

from optparse import (OptionParser,BadOptionError,AmbiguousOptionError)


class PassThroughOptionParser(OptionParser):

    def _process_args(self, largs, rargs, values):
        while rargs:
            try:
                self._process_args2(largs,rargs,values) 
            except (BadOptionError,AmbiguousOptionError), e:  # On failure, pass option to output list
                if sys.version_info < (2, 6, 0):
                    # Port to Python 2.4
                    p = re.match("^.*?no such option:\s*(.*?)$", e.msg)
                    if p:
                        largs.append(p.group(1))
                else:
                    largs.append(e.opt_str)
    

    # This version of the function successfully passes through negative numbers    
    def _process_args2(self, largs, rargs, values):
        """_process_args(largs : [string],
                         rargs : [string],
                         values : Values)

        Process command-line arguments and populate 'values', consuming
        options and arguments from 'rargs'.  If 'allow_interspersed_args' is
        false, stop at the first non-option argument.  If true, accumulate any
        interspersed non-option arguments in 'largs'.
        """
        while rargs:
            arg = rargs[0]

            p = re.match('^-[0-9.]+$', arg)
            if p:
                del rargs[0]
                raise BadOptionError(arg)
                #self.error(_("%s unrecognized number in arguments") % arg)
            
            # We handle bare "--" explicitly, and bare "-" is handled by the
            # standard arg handler since the short arg case ensures that the
            # len of the opt string is greater than 1.
            if arg == "--":
                del rargs[0]
                return
            elif arg[0:2] == "--":
                # process a single long option (possibly with value(s))
                OptionParser._process_long_opt(self, rargs, values)
            elif arg[:1] == "-" and len(arg) > 1:
                # process a cluster of short options (possibly with
                # value(s) for the last one only)
                OptionParser._process_short_opts(self, rargs, values)
            elif self.allow_interspersed_args:
                largs.append(arg)
                del rargs[0]
            else:
                return                  # stop now, leave this arg in rargs

        # Say this is the original argument list:
        # [arg0, arg1, ..., arg(i-1), arg(i), arg(i+1), ..., arg(N-1)]
        #                            ^
        # (we are about to process arg(i)).
        #
        # Then rargs is [arg(i), ..., arg(N-1)] and largs is a *subset* of
        # [arg0, ..., arg(i-1)] (any options and their arguments will have
        # been removed from largs).
        #
        # The while loop will usually consume 1 or more arguments per pass.
        # If it consumes 1 (eg. arg is an option that takes no arguments),
        # then after _process_arg() is done the situation is:
        #
        #   largs = subset of [arg0, ..., arg(i)]
        #   rargs = [arg(i+1), ..., arg(N-1)]
        #
        # If allow_interspersed_args is false, largs will always be
        # *empty* -- still a subset of [arg0, ..., arg(i-1)], but
        # not a very interesting subset!


