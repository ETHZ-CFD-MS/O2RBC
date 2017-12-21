#!/usr/bin/env python
#
# Extract the statistics of values in a postprocessing file produced by OpenFOAM
#

import os
import numpy as np
import argparse

from readPostProcessingFile import postProcessingFilePath, \
                                   readPostProcessingFile

def extractStatistics(postProcessingDict, fromTime = 0.0):

    times   = postProcessingDict['values'][:,0]
    values  = postProcessingDict['values'][:,1:]

    sampleTimeIdx = [i for i,t in enumerate(times) if t >= fromTime]

    means = np.mean(values[sampleTimeIdx,:], axis=0)

    print 'Means = %f' % means


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dirName', '-d', help='Name of the postprocessing directory')
    parser.add_argument('--fileName', '-f', help='File to load')
    parser.add_argument('--fromTime', type=float, help='Time from which extract statistics', 
                        default=0.0)

    args            = parser.parse_args()
    dirName         = args.dirName
    fileName        = args.fileName
    fromTime        = args.fromTime

    filePath = postProcessingFilePath('.', dirName, fileName)
    postProcessingDict = readPostProcessingFile(filePath)
    extractStatistics(postProcessingDict, fromTime)





