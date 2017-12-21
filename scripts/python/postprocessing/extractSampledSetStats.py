#!/usr/bin/env python
#
# Extract statistics from OpenFOAM sampled sets.

import numpy as np
import argparse

from postprocessing.loadSampledSetFiles import loadSampledSets


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

## Compute the time-averaged value of a sampled set.
#  @param sampledSetsDict Dictionary produced by loadSampledSets
def timeAveragedSampledSets(sampledSetsDict, min_time=0.0, max_time=1e100):

    positions = sampledSetsDict['positions']
    nPositions = np.size(positions, 0)
    nTimes = len(sampledSetsDict.keys()) - 1
    times = np.zeros((nTimes, 1))
    values = np.zeros((nPositions, nTimes))

    i = 0
    for key in sorted(sampledSetsDict):
        if is_number(key):
            time = float(key)
            times[i] = time
            values[:,i] = sampledSetsDict[key]
            i += 1

    sampleTimeIdx = [i for i,t in enumerate(times) if t >= min_time and t <= max_time]

    # currently assume a fixed time step
    averagedValues = np.mean(values[:,sampleTimeIdx], 1)

    return averagedValues


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sampledSet', '-s', help='Name of the sampledSet directory')
    parser.add_argument('--fileName', '-f', help='Name of the sampledSet files')
    parser.add_argument('--fromTime', type=float, help='Time from which extract statistics', 
                        default=0.0)

    args          = parser.parse_args()
    sampledSetDir = args.sampledSet
    fileName      = args.fileName
    fromTime     = args.fromTime

    sampledSetsDict = loadSampledSets('.', sampledSetDir, fileName)
    # sampledSetsDict = pickle.load(open('sampledSetDicts.pkl', 'r'))
    averagedValues = timeAveragedSampledSets(sampledSetsDict, fromTime)

    print 'Averaged values: '
    print averagedValues


