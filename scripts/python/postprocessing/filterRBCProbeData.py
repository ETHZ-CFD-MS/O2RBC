#!/usr/bin/env python
#
# Output the data recorded by a RBC probe with given edge index and 
# s-coordinate
#

import os
import numpy as np
import argparse

from loadRBCProbes import loadRBCProbes
from extractRBCProbeStats import RBCProbeValues, \
                                 listProbes 

separator = ', '

def outputRBCProbeData(probeDict, fromTime):
    headers = ['hitTimes', 'RBCIndex', 'meanValues', \
               'minValues', 'maxValues']
    formats = ['%g', '%i', '%g', '%g', '%g']
    print separator.join(headers)
    for i in range(len(probeDict['hitTimes'])):
        if probeDict['hitTimes'][i] >= fromTime:
            data = [probeDict[key][i] for key in headers]
            outputStrings = [f % x for f, x in zip(formats, data)]
            print separator.join(outputStrings)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--field', '-f', help='Field to sample', default='Hb')
    parser.add_argument('--probeName', '-p', help='Name of the probe directory')
    parser.add_argument('--edgeIndex', '-e', help='Edge index of the probe to sample')
    parser.add_argument('--sCoord', '-s', help='Curvilinear coordinate of the probe to sample')
    parser.add_argument('--fromTime', type=float, help='Time from which extract statistics', 
                        default=0.0)

    args       = parser.parse_args()
    fieldName  = args.field
    probeName  = args.probeName
    edgeIndex  = int(args.edgeIndex)
    sCoord     = float(args.sCoord)
    fromTime   = args.fromTime

    RBCProbesDict = loadRBCProbes('.', probeName, fieldName)
    probeDict = RBCProbeValues(RBCProbesDict, edgeIndex, sCoord)
    outputRBCProbeData(probeDict, fromTime)
