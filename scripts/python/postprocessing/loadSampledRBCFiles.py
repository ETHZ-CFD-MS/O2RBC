"""
Load results of the sampleRBCField function object into a dictionary
"""

import cPickle as pickle
import fnmatch
import os

import numpy as np

from readPostProcessingFile import postProcessingDirPath, \
                                   postProcessingFilePaths, \
                                   readPostProcessingFile
from utilities.containers import merge_list_from_dict_list


def loadSampledRBCs(case_name, sample_dir_name, cache_paths=True):
    """
    Load RBC fields produced by sampleRBCField into a dictionary.

    Args:
        case_name (str): case name
        sample_dir_name (str): directory name that contains the sampled fields
        cache_paths (bool): store the loaded RBCs in a pickle file for faster reloading

    Returns:
        Dictionary indexed by RBC names that contains the RBC fields
    """
    pickleName = 'sampledRBCDicts.pkl'

    pathToPickle = os.path.join(case_name, pickleName)
    if cache_paths and os.path.exists(pathToPickle):
        sampledRBCDict = pickle.load(open(pathToPickle, 'r'))
        print "Read pickle file %s" % pathToPickle

    else:
        print 'Reading sampled RBCs in {:s}'.format(case_name)
        sampledRBCDict = {}
        pathToPostProcessing = postProcessingDirPath(case_name, sample_dir_name)

        # extract all the RBC names
        RBCNames = []
        for root, subdirs, files in os.walk(pathToPostProcessing):
            for fileName in files:
                if fnmatch.fnmatch(fileName, 'RBC*'):
                    RBCNames.append(fileName)
        RBCNames = list(set(RBCNames))

        for RBCName in RBCNames:
            pathsToRBCs = postProcessingFilePaths(case_name, sample_dir_name, RBCName)
            sampledRBCDict[RBCName] = readSampledRBCDictFromPaths(pathsToRBCs)

        if cache_paths:
            pickle.dump(sampledRBCDict, open(pathToPickle, 'wb'))
            print "Wrote sampled RBCs to %s" % pathToPickle

    return sampledRBCDict


def readSampledRBCDictFromPaths(filePaths):
    """
    Return a dictionary with information for the sampled RBC from the corresponding
    file paths. The sampled RBC may be spread over several files.

    Args:
        filePaths (list): List of paths

    Returns:
        Dictionary with sampled RBC data

    """
    sampledRBCDicts = []
    for filePath in filePaths:
        sampledRBCDicts.append(readSampledRBCDictFromOneFile(filePath))
    return mergeSampledRBCDicts(sampledRBCDicts)


def mergeSampledRBCDicts(sampledRBCDicts):
    mergedDict = {key: merge_list_from_dict_list(key, sampledRBCDicts) for key in sampledRBCDicts[0]}
    unique_times, unique_ids = np.unique(mergedDict['time'], return_index=True)
    for key in mergedDict:
        mergedDict[key] = mergedDict[key][unique_ids]
    return mergedDict


def readSampledRBCDictFromOneFile(filePath):
    """
    Create a dict with data for the RBC data stored in filePath.

    This reader handles two different formats for sampled RBCs.
    One format is for simulations in graphs where the RBC position is given by a curvilinear
    coordinate with keys 'edgeIndex' and 'sCoord'.
    The other format is for simulations in axisymmetric domains where the longitudinal
    position is given by the key 'x'.

    Args:
        filePath (str): path to file

    Returns:
        dict with keys as below
    """
    postProcessingDict = readPostProcessingFile(filePath)
    if 'EdgeIndex' in postProcessingDict['header']:
        return {'time': postProcessingDict['values'][:,0],
                'mean': postProcessingDict['values'][:,1],
                'min': postProcessingDict['values'][:,2],
                'max': postProcessingDict['values'][:,3],
                'edgeIndex': postProcessingDict['values'][:,4].astype(int),
                'sCoord': postProcessingDict['values'][:,5]}
    else:
        return {'time': postProcessingDict['values'][:,0],
                'x': postProcessingDict['values'][:,1],
                'mean': postProcessingDict['values'][:,2],
                'min': postProcessingDict['values'][:,3],
                'max': postProcessingDict['values'][:,4]}
