#!/usr/bin/env python
#
# Load all x-profiles in a multi-RBC simulation.

import os
import numpy as np

## Load probes produced by OpenFOAM
#  @param caseName Path to current case
#  @param probeName Name of the probe directory
#  @param fieldName Name of the field to load
#  @return probes Dictionary with keys 'positions', 'nProbes', 'name', 'fieldName', 'times', 'values'.
def loadProbes(caseName, probeName, fieldName):

    pathToProbe = '%s/%s/0/%s' % (caseName, probeName, fieldName)
    probePositions = []
    # read the number of lines
    n_steps = sum(1 for line in open(pathToProbe)) - 4
    times = np.zeros((n_steps))
    # read the header
    with open(pathToProbe, 'r') as f:
        i = 0
        values = np.zeros((0,0))
        n_probes = 0

        for line in f:
            i = i+1
            if i == 1:
                n_probes = len(line.split()) - 2
                probePositions = np.zeros((3, n_probes))
            if i <= 3:
                coords = line.split()[2:]
                probePositions[i-1,:] = coords
            elif i == 4:
                values = np.zeros((n_steps, n_probes))
            elif i >= 5:
                t = i-5
                split_line = [float(x) for x in line.split()]
                times[t] = split_line[0]
                values[t,:] = split_line[1:]

    return dict([('positions', probePositions), \
                 ('name', probeName), \
                 ('fieldName', fieldName), \
                 ('nProbes', n_probes), \
                 ('times', times), \
                 ('values', values)])

if __name__ == "__main__":
    loadProbes('domain', 'probeMidstreamPO2', 'PO2')

