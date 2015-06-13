#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, glob, optparse, re, shutil, subprocess, string, time, logging, threading


def callGrep(name, location):

    # Call grep to get a list of possible locations
    cmd = ['grep', '-n', '-r', '-I', '--include', '*.h', '--include', '*.hpp', '--include', '*.tcc', '-w', name, location]
    #print cmd
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    cmdText, err = p.communicate()

    # Now filter out lines which are probably not a definition    
    outputText = ''
    for line in cmdText.split('\n'):
        
        if len(line) < len(name): # Skip obvious non-matches
            continue
        
        codeStart = line.find(':')
        lineMatch = line.find(name, codeStart+1)
        
        # If we see any of these before the name on this line it probably is not a definition
        skipThis = False
        banList = ['=', '.', 'return', '//', 'std::cout', 'if', '->']
        for b in banList:
            match = line.find(b, codeStart+1)
            if (match >= 0) and (match < lineMatch):
                skipThis = True
                break
        if skipThis:
            continue
        
        # If we see any of these before the name on this line it probably is a definition
        markThis = False
        markList = ['struct', 'class']
        for m in markList:
            if m in line:
                markThis = True
        if markThis: # Do something to make the line stand out a little
            outputText = outputText + '\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n'
        
        #print line
        #print line[:codeStart]
        #print line[lineMatch:]
        #print '---------'
        # Cut out the location text that was passed in to declutter the output
        outputText = outputText + line[len(location):] + '\n'
    
    return outputText

def findAspThing(name):
    
    # Search VW
    outputText = callGrep(name, '/home/smcmich1/repo/visionworkbench/src/vw')

    print '\n############# VW ##################################################'    
    print outputText
    
    # Search ASP
    outputText = callGrep(name, '/home/smcmich1/repo/StereoPipeline/src/asp')

    print '\n############# ASP ##################################################'    
    print outputText


def main(argsIn):

    findAspThing(argsIn[0])
    
    
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))