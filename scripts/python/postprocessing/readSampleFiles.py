"""Utilities to read sample files generated by OpenFOAM."""

import numpy as np
import sys
import os


def loadSampleFile(caseName, setDirectory, time, setName, fieldNames=['PO2']):
    """Return a numpy array with the sample file values."""
    pathToSample = createSampleFilePath(caseName, setDirectory, time,
                                        setName, fieldNames)
    num_lines = sum(1 for line in open(pathToSample))

    if num_lines==0:
        sys.exit("error: file %s has zero lines" % pathToSample)

    with open(pathToSample, 'r') as f:
        array = np.zeros( (num_lines, len(fieldNames)+1) )
        i = 0
        for line in f:
            array[i,:] = [float(x) for x in line.split()]
            i = i+1

    return array

def loadSampleFileFromTimes(caseName, setDirectory, times, 
                            setName, fieldNames=['PO2']):
    """Return a dictionary with keys 'x' and the elements of fieldNames.

    The values corresponding to the fieldNames are numpy arrays with shape
    (len(coords), len(times))
    """
    sampleDataStart = loadSampleFile(caseName, setDirectory,
                                     times[0], setName, fieldNames)
    coords = sampleDataStart[:,0]
    resultDict = {'x': coords}
    for fieldName in fieldNames:
        resultDict[fieldName] = np.zeros((len(coords), len(times)))

    for i, time in enumerate(times):
        sampleData= loadSampleFile(caseName, setDirectory,
                                   time, setName, fieldNames)
        for j, fieldName in enumerate(fieldNames):
            resultDict[fieldName][:,i] = sampleData[:,j+1]
    return resultDict

def createSampleFilePath(caseName, setDirectory, time, setName, fieldNames=['PO2']):
    setExtension = 'xy'
    timeString = '%g' % time
    fieldNameString = '_'.join(fieldNames)
    pathToSample = os.path.join(caseName, 'domain', 'postProcessing',
                                setDirectory, timeString,
                                '%s_%s.%s' % (setName, fieldNameString, setExtension))
    return pathToSample
