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

import sys

import os, glob, optparse, re, shutil, subprocess, string, time

#import git

import IsisTools, IrgFileFunctions, IrgIsisFunctions, IrgStringFunctions

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''Tool for managing lronac pipeline production data'''

    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg



#--------------------------------------------------------------------------------

# All status information is saved in this folder
VW_FOLDER      = '/home/smcmich1/repo/visionworkbench'
ASP_FOLDER     = '/home/smcmich1/repo/StereoPipeline'
TOOLS_FOLDER   = '/home/smcmich1/repo/Tools'
TOOLS_B_FOLDER = '/home/smcmich1/repo/ToolsBuild'
LRO_FOLDER     = '/home/smcmich1/repo/lronacPipeline'
LRO_B_FOLDER   = '/home/smcmich1/repo/lronacPipelineBuild'

#
#def updateRepo(repoFolder, buildFolder):
#    """Updates the contents of a repository"""
#    
#    # Account for in-place builds
#    if not buildFolder:
#        buildFolder = repoFolder
#    
#    # Set up repo object
#    repoName = os.path.basename(repoFolder)
#    repo     = git.Repo(repoFolder)
#
#    # Up to the user to handle local changes    
#    if repo.is_dirty():
#        raise Exception('Cant operate on dirty repo ' + repoName)
#    
#    
#    # Find the upstream and origin remotes and fetch
#    upstreamRemote = None
#    originRemote   = None
#    for r in repo.remotes:
#        if r.name == 'upstream':
#            upstreamRemote = r
#        if r.name == 'origin':
#            originRemote = r
#
#    if not originRemote: # Upstream remote is optional
#        raise Exception('Failed to find origin remote!')
#
#    if upstreamRemote:
#        upstreamRemote.pull('--rebase')
#        #fetchResult = r.fetch()
#        #for f in fetchResult:
#        #    print f
#            
#            
#    

class GitRepo:
    """Class for working with a git repository"""
    
    _gitFolder  = ''
    _repoFolder = ''
    
    def __init__(self, folder):
        """Constructor"""
        self._repoFolder = folder
        self._gitFolder  = os.path.join(folder, '.git/')
        
    def callGitCommand(self, gitCommandList):
        """Calls a set of git commands on a folder and returns the output text"""
    
        # Get the tag using a subprocess call
        cmd = ['git', '--git-dir', self._gitFolder, '--work-tree', self._repoFolder]
        cmd = cmd + gitCommandList
        #print cmd
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        textOutput, err = p.communicate()
        
        return textOutput
    
    
    def getBranch(self):
        """Returns the name of the current branch"""
        textOutput = self.callGitCommand(['branch'])
        for line in textOutput.split('\n'):
            if line.startswith('*'):
                return line[2:]
        raise Exception('Unable to determine branch!')
    
    def isDirty(self):
        """Returns True if the repo is dirty"""    
        
        textOutput = self.callGitCommand(['status'])
        return (textOutput.find('Changes not staged for commit') >= 0)

    
    def hasUpstreamRemote(self):
        """Returns True the repo has an upstream source"""    
        
        textOutput = self.callGitCommand(['remote'])   
        return (textOutput.find('upstream') >= 0)
    

    
    
def updateRepo(repoFolder, buildFolder):
    """Updates the contents of a repository"""
    
    print 'Updating repository: ' + repoFolder
    
    # Account for in-place builds
    if not buildFolder:
        buildFolder = repoFolder
    
    # Set up repo object
    repo = GitRepo(repoFolder)
        
    if repo.isDirty():
        raise Exception('TODO: Repo dirty! Add stash functionality!')
    
    if repo.getBranch() != 'master':
        raise Exception('TODO: Not master branch! Add branch functionality!')
    
    fetchedAnything = False
    if repo.hasUpstreamRemote():
        print 'Fetching updates from upstream...'
        # Fetch from the upstream remote if it exists
        upstreamFetchText = repo.callGitCommand(['fetch', 'upstream'])
        
        # If we pulled down anything new, rebase on to it
        if len(upstreamFetchText) > 10: # If we fetched anything
            fetchedAnything = True
            print 'Rebasing on to upstream updates...'
            upstreamRebaseText = repo.callGitCommand(['rebase', 'upstream/master'])
            if (upstreamRebaseText.find('conflict') >= 0):
                print upstreamRebaseText
                raise Exception('Found merge conflict in upstream rebase!')
    
   
    # Fetch from the origin remote
    print 'Fetching updates from origin...'
    originFetchText = repo.callGitCommand(['fetch', 'origin'])
    
    # If we pulled down anything new, rebase on to it
    if len(originFetchText) > 10: # If we fetched anything
        fetchedAnything = True
        print 'Rebasing on to origin updates...'
        originRebaseText = repo.callGitCommand(['rebase', 'origin/master'])
        if (originRebaseText.find('conflict') >= 0):
                print originRebaseText
                raise Exception('Found merge conflict in origin rebase!')


    # If we didn't get any new code we are finished!
    if not fetchedAnything:
        print 'Did not fetch any new updates, finished with ' + repoFolder
        return True
    
    # At this point our code should be up to date, try to build
    # Move to build folder
    cmd = 'cd ' + buildFolder
    os.system(cmd)

    # Build!
    cmd = ['make', '-j', '8']
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    makeOutput, err = p.communicate()
    if (makeOutput.find('make: *** [all] Error') >= 0):
        raise Exception('Error building on folder ' + buildFolder)
    
    # TODO: Don't do this if we didn't get new code
    if buildFolder == repoFolder: # This is a hack for my repositories!
        # Try calling installation
        cmd = 'make install '
        os.system(cmd)
        # TODO: Is there a way to check this?
    
    
    
    return True
    
    
    

def main():

    #
    #try:
    #    try:
    #        usage = "usage: archiveManager.py [--help][--manual]\n  "
    #        parser = optparse.OptionParser(usage=usage)
    #        
    #        parser.add_option("--find-pbs", dest="pbsName",
    #                          help="Get the full path to the data set for a PBS job.",
    #                          default=None)
    #
    #
    #        parser.add_option("--flag-incomplete", action="store_true",
    #                          dest="flagIncomplete", default=False,
    #                          help="Flag incomplete data sets for repeat processing.")
    #        
    #        parser.add_option("--clear", action="store_true", dest="clear", default=False,
    #                          help="Clear files after archiving them.")
    #        
    #        parser.add_option("--dry-run", action="store_true", dest="dryRun", default=False,
    #                          help="Don't touch any data, just display actions to be taken.")
    #                          
    #        parser.add_option("--manual", action="callback", callback=man,
    #                          help="Read the manual.")
    #        (options, args) = parser.parse_args()
    #
    #    except optparse.OptionError, msg:
    #        raise Usage(msg)


    repoList = [(VW_FOLDER,  None),
                (ASP_FOLDER, None),
                (TOOLS_FOLDER, TOOLS_B_FOLDER),
                (LRO_FOLDER, LRO_B_FOLDER)]
    
    for r in repoList:
        updateRepo(r[0], r[1])
    


    return 0

    #except Usage, err:
    #    print >>sys.stderr, err.msg
    #    return 2


if __name__ == "__main__":
    sys.exit(main())
