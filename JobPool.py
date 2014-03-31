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

"""JobPool.py - Simple functions for managing parallel threaded tasks"""

import sys, os, re, shutil, subprocess, string, time, errno

# This module could probably be improved a lot!

job_pool  = [] # currently running jobs


def add_job(cmd, numProcesses, verbose=True):
    """Adds a job to the pool when space is available"""
    SLEEP_TIME = 0.001
    
    while ( len(job_pool) >= numProcesses ): # Wait until another process finishes
        for i in range(len(job_pool)):
            if ( job_pool[i].poll() is not None ):
                job_pool.pop(i)
                job_pool.append( subprocess.Popen(cmd) )
                __add_job_internal(cmd, verbose)
                return
        time.sleep( SLEEP_TIME )
        sleep_time = (SLEEP_TIME * 5) % 60
        
    __add_job_internal(cmd, verbose) # Did not have to wait, just add the job


def __add_job_internal(cmd, verbose):
    """Add a job to the pool without regard to room"""
    
    if not verbose:
        FNULL = open(os.devnull, 'w')
        job_pool.append( subprocess.Popen(cmd, stdout=FNULL, stderr=subprocess.STDOUT) )
    else: # verbose
        job_pool.append( subprocess.Popen(cmd) )


def wait_on_all_jobs():
    """Wait until all queued jobs are finished before exiting this function"""
    
    SLEEP_TIME = 1
    while len(job_pool) > 0:
        for i in range(len(job_pool)):
            if ( job_pool[i].poll() is not None ):
                job_pool.pop(i)
                break # must restart as array changed size
        time.sleep( SLEEP_TIME )




