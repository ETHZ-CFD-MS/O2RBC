#!/usr/bin/env python
#
# Extract the statistics of values from the function object RBCProbe
#

import os
import numpy as np
import argparse

from loadRBCProbes import loadRBCProbes


## Extract the RBC probe statistics at a given graph position
#
# @param caseName Name of the case where probe data are stored
# @param RBCProbeDict Dictionary produces by loadRBCProbes
# @param edgeIndex Edge index
# @param sCoord Curvilinear coordinate
# @return Dictionary with keys 'meanValues', 'minValues' and 'maxValues'
def RBCProbeValues(RBCProbesDict, edgeIndex, sCoord):
    returnDict = {'times':      [], \
                  'hitTimes':   [], \
                  'RBCIndex':   [], \
                  'meanValues': [], \
                  'minValues':  [], \
                  'maxValues':  []}

    for i in range(len(RBCProbesDict['times'])):
        if     RBCProbesDict['edgeIndices'][i] == edgeIndex \
          and (RBCProbesDict['sCoords'][i]     == sCoord):
            for key in returnDict:
                returnDict[key].append(RBCProbesDict[key][i])

    for key in returnDict:
        returnDict[key] = np.asarray(returnDict[key])

    if len(returnDict['times']) == 0:
        raise ValueError('No entry found for edgeIndex %i and '
                         'sCoord %g' % (edgeIndex, sCoord))

    return returnDict


## Extract statistics of sampled RBC values for all probes
def extractAllRBCProbeStatistics(RBCProbesDict, fromTime=0.0):
    probeList = listProbes(RBCProbesDict)
    statDicts = []
    for (eI, s) in probeList:
        statDicts.append(extractRBCProbeStatistics(RBCProbesDict, eI, s, fromTime))
    headers = ['edgeIndex', 'sCoord']
    headers.extend(sortedStatDictKeys(statDicts[0]))
    print ', '.join(headers)
    for i, (eI, s) in enumerate(probeList):
        stats = [statDicts[i][key] for key in sortedStatDictKeys(statDicts[i])]
        outputStrings = ['%i' % eI, '%g' % s]
        outputStrings.extend(['%g' % x for x in stats])
        print ', '.join(outputStrings)
        

## Extract statistics of sampled RBC values for a given probe
def extractRBCProbeStatistics(RBCProbesDict, edgeIndex, sCoord, fromTime=0.0):

    probeDict = RBCProbeValues(RBCProbesDict, edgeIndex, sCoord)

    times      = probeDict['times']
    meanValues = probeDict['meanValues']
    minValues  = probeDict['minValues']
    maxValues  = probeDict['maxValues']

    sampleTimeIdx = [i for i,t in enumerate(times) if t >= fromTime]

    meanMean = [np.mean(meanValues[sampleTimeIdx])]
    meanMin  = [np.mean(minValues[sampleTimeIdx])]
    meanMax  = [np.mean(maxValues[sampleTimeIdx])]
    stdMean = [np.std(meanValues[sampleTimeIdx])]
    stdMin  = [np.std(minValues[sampleTimeIdx])]
    stdMax  = [np.std(maxValues[sampleTimeIdx])]

    return_dict = \
            dict([('meanMean', meanMean), \
                  ('meanMin',  meanMin), \
                  ('meanMax',  meanMax), \
                  ('stdMean', stdMean), \
                  ('stdMin',  stdMin), \
                  ('stdMax',  stdMax)])

    # convert to numpy arrays
    for key in return_dict.keys():
        return_dict[key] = np.asarray(return_dict[key])

    return return_dict

def sortedStatDictKeys(statDict):
    return sorted(statDict.keys())

## Return the list of all probes contained in an RBCProbesDict
#  @return List of tuples (edgeIndex, sCoord)
def listProbes(RBCProbesDict):
    tupleList = [(eI, s) for eI, s in zip(RBCProbesDict['edgeIndices'], RBCProbesDict['sCoords']    )]
    return sorted(list(set(tupleList)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--field', '-f', help='Field to sample', default='Hb')
    parser.add_argument('--probeName', '-p', help='Name of the probe directory')
    # parser.add_argument('--edgeIndex', '-e', help='Edge index of the probe to sample')
    # parser.add_argument('--sCoord', '-s', help='Curvilinear coordinate of the probe to sample')
    parser.add_argument('--fromTime', type=float, help='Time from which extract statistics', 
                        default=0.0)

    args            = parser.parse_args()
    fieldName       = args.field
    probeName       = args.probeName
    # edgeIndex       = int(args.edgeIndex)
    # sCoord          = float(args.sCoord)
    fromTime       = args.fromTime

    RBCProbesDict = loadRBCProbes('.', probeName, fieldName)
    extractAllRBCProbeStatistics(RBCProbesDict, fromTime)

    # RBCStatDict = extractRBCProbeStatistics(RBCProbesDict, edgeIndex, sCoord, fromTime)

    # positionStr = '(%i, %5.2e)' % (edgeIndex, sCoord)
    # print 'RBC probe at position %s:' % positionStr
    # for key, value in sorted(RBCStatDict.iteritems()):
        # print '%s = %g' % (key, value)

