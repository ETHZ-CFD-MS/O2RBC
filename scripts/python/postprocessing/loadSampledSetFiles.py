"""
Load results from OpenFOAM sampled sets into a dictionary
"""

import os
import numpy as np
import cPickle as pickle

from readPostProcessingFile import postProcessingDirPath, \
                                   postProcessingTimeSteps, \
                                   readPostProcessingFile

def loadSampledSets(caseName, sampledDirName, sampledSetFile,
        maxTimeInterval=float('inf')):
    """
    Load sampled sets produced by OpenFOAM. Caches results in a pickle file.

    If maxTimeInterval is smaller than the time span covered by the sampled files,
    only the time steps between lastTime - maxTimeInterval and lastTime are loaded.

    Args:
        caseName: Path to current case
        sampledDirName: Name of the sample sets directory
        sampledSetFile: Name of the sample set files
        maxTimeInterval: Maximal time interval which is loaded.
    Returns:
        Dictionary with keys corresponding to the time values
        of the output time steps, and one key 'position'
        for the point positions
    """
    pickleName = 'sampledSetDicts.pkl'
    pathToPickle = os.path.join(caseName, pickleName)
    try:
        sampledSetsDict = pickle.load(open(pathToPickle, 'r'))
        print "Read pickle file %s" % pathToPickle
    except IOError:
        sampledSetsDict = {}
        sampleDirPath = postProcessingDirPath(caseName, sampledDirName)
        print "Reading sampled sets..."
        timeDirs = postProcessingTimeSteps(caseName, sampledDirName)
        lastTime = max([float(t) for t in timeDirs])
        timeDirs = filter(lambda t: lastTime - float(t) <= maxTimeInterval, timeDirs)
        for timeDir in timeDirs:
            filePath = os.path.join(sampleDirPath, timeDir, sampledSetFile)
            fullDict = sampledSetDict(filePath)
            sampledSetsDict[timeDir] = fullDict['values']

        sampledSetsDict['positions'] = fullDict['positions']
        pickle.dump(sampledSetsDict, open(pathToPickle, 'wb'))
        print "Wrote sampled sets to %s" % pathToPickle
    return sampledSetsDict

def sampledSetDict(filePath):
    """
    Return a dictionary with sampledSet information from the file path
    of that sampledSet
    """
    postProcessingDict = readPostProcessingFile(filePath)
    try:
        return {'positions': postProcessingDict['values'][:,0],
                'values':    postProcessingDict['values'][:,1]}
    except IndexError:
        print 'Index error while loading file {}'.format(filePath)
        raise

