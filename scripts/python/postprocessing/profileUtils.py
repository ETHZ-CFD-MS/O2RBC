#!/usr/bin/env python
#
# Utilities from OpenFOAM probes.
#

import numpy as np

from parse.readRBCPositions import readRBCPositions


## Return times at which some key part of a RBC (front, hollow or rear)
#  crosses the x-coordinate of a probe
# @param time Current time
# @param RBCPositionsPath Path to file with RBC positions
# @param RBCPart Part of the RBC: 'front', 'hollow' or 'rear'
# @return Numpy array with corresponding positions
def getRBCPositions(time, RBCPositionsPath, RBCPart):

    # load RBC positions
    (RBCPositions, times) = readRBCPositions(RBCPositionsPath)
    n_RBC = len(RBCPositions)

    # find time index
    t_idx = np.searchsorted(times, time)

    # for each RBC, find position at that time
    positions = []
    for i in range(n_RBC):
        RBCFronts = RBCPositions[i][RBCPart]
        positions.append(RBCFronts[t_idx])

    positions.sort()
    positions = np.asarray(positions)
    return positions
    







