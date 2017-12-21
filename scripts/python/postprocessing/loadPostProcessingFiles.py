#!/usr/bin/env python
#
# Load results of OpenFOAM probes into a dictionary
# Utilities to read postprocessing files

import os
import sys
import warnings

import numpy as np

from utilities.sort import natural_keys

sampleDir = 'sets'
postProcessingDir = 'postProcessing'

## Find paths to postprocessing files and returns them in a list
#  @param caseName Path to the current case
#  @param dirName Name of the postprocessing directory
#  @param fileName Name of the files to load
#  @return pathsToProbes list with the paths to the postprocessing files
def getPathsToPostProcessingFiles(caseName, dirName, fileName):
    if os.path.isdir('%s/postProcessing' % caseName):
        dirList = readDirectories('%s/postProcessing/%s' % (caseName, dirName))
        checkPostProcessingDirectories(dirList)
        if len(dirList) > 1:
            print 'Loading postprocessing files from directories %s' % ', '.join(dirList)
        pathsToProbes = ['%s/postProcessing/%s/%s/%s' % (caseName, dirName, dir, fileName) for dir in dirList]
    else:
        dirList = readDirectories('%s/%s' % (caseName, dirName))
        checkPostProcessingDirectories(dirList)
        if len(dirList) > 1:
            print 'Loading postprocessing files from directories %s' % ', '.join(dirList)
        pathsToProbes = ['%s/%s/%s/%s' % (caseName, dirName, dir, fileName) for dir in dirList]

    pathsToProbes = [path for path in pathsToProbes if os.path.isfile(path)]
    if pathsToProbes == []:
        message = 'WARNING, found no file for caseName = %s, dirName = %s and fileName = %s' \
                        % (caseName, dirName, fileName)
        warnings.warn(message)

    return pathsToProbes

## Create path to a file produced by the OpenFOAM utility 'sample'
def createSampleFilePath(caseName, time, sampleFile):
    pathToSample = os.path.join(caseName, postProcessingDir, sampleDir, time, sampleFile)
    return pathToSample

## Load a sample file produced by the OpenFOAM utility 'sample'
def loadSampleFile(caseName, time, sampleFile):
    pathToSample = createSampleFilePath(caseName, time, sampleFile)
    num_lines = sum(1 for line in open(pathToSample))
    if num_lines==0:
        sys.exit("error: file %s has zero lines" % pathToSample)

    with open(pathToSample, 'r') as f:
        i = 0
        for line in f:
            if i == 0:
                nCol = len(line.split())
                array = np.zeros( (num_lines, nCol ) )

            array[i,:] = [float(x) for x in line.split()]
            i = i+1

    return array

def readDirectories(path):
    dirList = os.walk(path).next()[1]
    dirList.sort(key=natural_keys)
    return dirList

def checkPostProcessingDirectories(dirNames):
    for d in dirNames:
        if not is_number(d):
            raise ValueError('The directory %s in the postprocessing directory is not numeric' % d)

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

