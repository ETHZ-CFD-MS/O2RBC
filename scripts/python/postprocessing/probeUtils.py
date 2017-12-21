#!/usr/bin/env python
#
# Utilities from OpenFOAM probes.
#

import os

from parse.readRBCPositions import readRBCPositions


# return times at which some key part of a RBC (front, hollow or rear)
# crosses the x-coordinate of a probe
def getRBCPassingTimes(probePosition, RBCPositionsPath, RBCPart):
    xProbe = probePosition[0]
    (RBCPositions, times) = readRBCPositions(RBCPositionsPath)
    n_RBC = len(RBCPositions)
    # for each RBC, find out crossing times
    cross_times = []
    for i in range(n_RBC):
        RBCFronts = RBCPositions[i][RBCPart]
        for j in range(len(times)-1):
            if xProbe >= RBCFronts[j] and xProbe < RBCFronts[j+1]:
                # find crossing time exactly using a linear
                # interpolation
                t = times[j] + (times[j+1] - times[j]) \
                             * (xProbe - RBCFronts[j]) \
                             / (RBCFronts[j+1] - RBCFronts[j])
                cross_times.append(t)

    cross_times.sort()
    return cross_times

# get index of closest RBC
def getClosestRBC(probePosition, RBCPositions, t_index):

    xProbe = probePosition[0]
    n_RBC = len(RBCPositions)

    # for each RBC, find distance to probe
    RBCCenters = [RBCPositions[i]['center'][t_index] for i in range(n_RBC)]
    i_closest = min(range(n_RBC), key=lambda i: abs(RBCCenters[i] - xProbe))
    # print i_closest, abs(RBCCenters[i_closest] - xProbe)

    return i_closest

## Return the names of all directories that start with 'probe' in the given folder
#  @param path Path to search in
#  @return probes List with names of probes
def probeNames(path, suffix=''):
    probes = []
    for filename in os.listdir(path):
        if filename.startswith('probe') and os.path.isdir('%s/%s' % (path, filename)) \
           and filename.endswith(suffix):
            probes.append(filename)

    probes.sort()

    return probes


if __name__ == "__main__":
    getRBCPassingTimes([60e-6, 0, 0], 'RBCPositions.txt', 'front')



    






