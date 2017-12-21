"""
Load results of the sampleRBCField function object into a dictionary
"""

import os

import numpy as np

from loadPostProcessingFiles import getPathsToPostProcessingFiles
from utilities.containers import check_equal, merge_list_from_dict_list


def loadRBCProbes(caseName, probeName, fieldName):
    """
    Args:
        caseName (str): Path to current case
        probeName (str): Name of the probe directory
        fieldName (str): Name of the field to load

    Returns:
        Dictionary with keys 'times', 'fieldName', 'edgeIndices', 'sCoords',
        'meanValues', 'minValues', 'maxValues'
    """
    pathsToProbes = getPathsToPostProcessingFiles(os.path.join(caseName, 'domain'),
                                                  probeName, fieldName)
    RBCProbeDicts = []
    for pathToProbe in pathsToProbes:
        RBCProbeDicts.append(loadRBCProbeFromOneFile(pathToProbe, fieldName))
    return mergeRBCProbeDicts(RBCProbeDicts)


def loadRBCProbeFromOneFile(pathToProbe, fieldName):
    n_steps = sum(1 for line in open(pathToProbe)) - 1
    with open(pathToProbe, 'r') as f:
        i = 0
        times = np.zeros((n_steps))
        edgeIndices = np.zeros((n_steps))
        sCoords     = np.zeros((n_steps))
        RBCIndex    = np.zeros((n_steps))
        hitTimes    = np.zeros((n_steps))
        meanValues  = np.zeros((n_steps))
        minValues   = np.zeros((n_steps))
        maxValues   = np.zeros((n_steps))

        for line in f:
            i = i+1
            if i >= 2:
                t = i-2
                split_line = [float(x) for x in line.split()]
                times[t]       = split_line[0]
                edgeIndices[t] = split_line[1]
                sCoords[t]     = split_line[2]
                RBCIndex[t]    = split_line[3]
                hitTimes[t]    = split_line[4]
                meanValues[t]  = split_line[5]
                minValues[t]   = split_line[6]
                maxValues[t]   = split_line[7]

    return {'times': times,
            'fieldName': fieldName,
            'edgeIndices': edgeIndices,
            'sCoords': sCoords,
            'hitTimes': hitTimes,
            'RBCIndex': RBCIndex,
            'meanValues': meanValues,
            'minValues': minValues,
            'maxValues': maxValues}


def mergeRBCProbeDicts(probeDicts):
    check_equal([probeDict['fieldName'] for probeDict in probeDicts])
    fieldName = probeDicts[0]['fieldName']
    times = np.append(probeDicts[0]['times'],
                      [probeDict['times'] for probeDict in probeDicts[1:]])
    returnDict = {'times': times,
                  'fieldName': fieldName}
    for key in probeDicts[0]:
        if key not in ['times', 'fieldName']:
            returnDict[key] = merge_list_from_dict_list(key, probeDicts)
    return returnDict

