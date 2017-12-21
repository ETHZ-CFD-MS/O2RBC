#!/usr/bin/env python
"""
Plot transit times and hemoglobin saturation of sampled RBCs.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from loadSampledRBCFiles import loadSampledRBCs
from plot.figureoptions import FigureOptions

HbInit = 0.73166

# style maps for 1st CTH graph
# traveledEdgeIDs = [68, 69, 70]
# traveledEdgeColors = ['r', 'g', 'b']
# edgeColorMap = dict([(68, 'r'), \
                     # (69, 'g'), \
                     # (70, 'b')])
# edgeStyleMap = dict([(65, '-'), \
                     # (66, '--'), \
                     # (67, ':')])
# markerStyleMap = dict([(65, 'o'), \
                       # (66, 's'), \
                       # (67, '^')])

# pooling color map for loop topology
traveledEdgeIDs = [3, 4]
traveledEdgeColors = ['r', 'b']
edgeStyleMap = dict([(3, '-'), \
                     (4, ':')])
edgeColorMap = dict([(3, 'r'), \
                     (4, 'b')])
markerStyleMap = dict([(3, 'o'), \
                       (4, 's')])

def sampledRBCLineStyle(sampledValues):
    edges = sampledValues['edgeIndex']

    for edgeIndex, style in edgeStyleMap.iteritems():
        if edgeIndex in edges:
            return style
    
    return '-'

def sampledRBCMarkerStyle(sampledValues):
    edges = sampledValues['edgeIndex']
    for edgeIndex, style in markerStyleMap.iteritems():
        if edgeIndex in edges:
            return style
    
    return 'o'

def sampledRBCColor(sampledValues):
    edges = sampledValues['edgeIndex']
    for edgeIndex, color in edgeColorMap.iteritems():
        if edgeIndex in edges:
            return color
    
    return 'k'

# Returns a dictionary with sampled RBCs that start after min_time and finish before max_time
def filterRBCs(sampledRBCsDict,  min_time, max_time):
    filteredDict = dict()
    for RBCName, sampledValues in sampledRBCsDict.iteritems():
        initialTime = sampledValues['time'][0]
        finalTime   = sampledValues['time'][-1]
        if initialTime >= min_time and finalTime <= max_time:
            filteredDict[RBCName] = sampledValues

    return filteredDict

def extractDataFromSampledRBCs(filteredDict):
    transitTimes = []
    finalHb = []
    deltaHb = []
    for RBCName, sampledValues in filteredDict.iteritems():
        transitTimes.append(sampledValues['time'][-1] - sampledValues['time'][0])
        finalHb.append(sampledValues['mean'][-1])
        deltaHb.append(HbInit - sampledValues['mean'][-1])

    return np.asarray(transitTimes), np.asarray(finalHb), np.asarray(deltaHb)

def extractDataPooledFromTraveledVessel(filteredDict):
    nPools = len(traveledEdgeIDs)
    transitTimes = [list() for _ in xrange(nPools)]
    finalHb = [list() for _ in xrange(nPools)]
    deltaHb = [list() for _ in xrange(nPools)]
    for RBCName, sampledValues in filteredDict.iteritems():
        for i, edgeI in enumerate(traveledEdgeIDs):
            if edgeI in sampledValues['edgeIndex']:
                transitTimes[i].append(sampledValues['time'][-1] - sampledValues['time'][0])
                finalHb[i].append(sampledValues['mean'][-1])
                deltaHb[i].append(HbInit - sampledValues['mean'][-1])

    return transitTimes, finalHb, deltaHb


def displayMoments(sampledRBCsDict, min_time, max_time):
    filteredDict = filterRBCs(sampledRBCsDict, min_time, max_time)
    transitTimes, finalHb, deltaHb = extractDataFromSampledRBCs(filteredDict)
    pooledTransitTimes, pooledFinalHb, pooledDeltaHb = extractDataPooledFromTraveledVessel(filteredDict)

    meanTT = np.mean(transitTimes)
    stdTT  = np.std(transitTimes)
    skewTT = stats.skew(transitTimes)
    rdTT   = stdTT/meanTT

    meanDeltaHb = HbInit - np.mean(finalHb)
    stdDeltaHb  = np.std(finalHb)
    skewDeltaHb = stats.skew(finalHb)
    rdDeltaHb   = stdDeltaHb/meanDeltaHb

    print 'Mean transit time = %g' % meanTT
    print 'Std  transit time = %g' % stdTT
    print 'Skew transit time = %g' % skewTT
    print 'Relative dispersion = %g' % rdTT

    print 'Mean final Hb = %g' % np.mean(finalHb)
    print 'Mean Delta Hb = %g' % meanDeltaHb
    print 'Std  final Hb = %g' % stdDeltaHb
    print 'Skew final Hb = %g' % skewDeltaHb
    print 'Relative dispersion of Delta Hb = %g' % rdDeltaHb

    nPools = len(pooledTransitTimes)
    pooledTT    = [np.mean(pooledTransitTimes[i]) for i in range(nPools)]
    pooledMeans = [np.mean(pooledFinalHb[i]) for i in range(nPools)]
    print
    print 'Pooled statistics:'
    print 'Mean transit time = ', pooledTT
    print 'Mean final Hb     = ', pooledMeans


def plotSampledRBCs(sampledRBCsDict, min_time, max_time):

    filteredDict = filterRBCs(sampledRBCsDict, min_time, max_time)
    transitTimes, finalHb, deltaHb = extractDataFromSampledRBCs(filteredDict)
    pooledTransitTimes, pooledFinalHb, pooledDeltaHb = extractDataPooledFromTraveledVessel(filteredDict)

    f, axarr = plt.subplots(2,2)#, sharex=True)

    # Plot [HbO] as function of transit time
    for RBCName, sampledValues in filteredDict.iteritems():
        initialTime = sampledValues['time'][0]
        finalTime   = sampledValues['time'][-1]
        transitTime = finalTime - initialTime

        color       = sampledRBCColor(sampledValues)
        linestyle   = sampledRBCLineStyle(sampledValues)
        markerstyle = sampledRBCMarkerStyle(sampledValues)

        plt.sca(axarr[0,0])
        plt.plot(sampledValues['time'] - initialTime, sampledValues['mean'], \
                 color=color, linestyle=linestyle)
        plt.plot(finalTime - initialTime, sampledValues['mean'][-1], \
                 color=color, marker=markerstyle)
        # plt.plot((sampledValues['time'] - initialTime)/transitTime, sampledValues['mean'], \
                 # color=color, linestyle=linestyle)
        # plt.plot((finalTime - initialTime)/transitTime, sampledValues['mean'][-1], \
                 # color=color, marker=markerstyle)
        plt.sca(axarr[1,0])
        plt.scatter(finalTime - initialTime, sampledValues['mean'][-1], \
                 edgecolor=color, facecolor='None', marker=markerstyle)

        if finalTime - initialTime < 0.01:
            print "RBC with short travel time: %s, %g, %g" % (RBCName, initialTime, finalTime)

    plt.sca(axarr[0,0])
    plt.xlabel(r'$\mathrm{Transit\; time}\; [s]$')
    plt.ylabel(r'$[\mathrm{HbO}]$')
    plt.sca(axarr[1,0])
    plt.xlabel(r'$\mathrm{Transit\; time}\; [s]$')
    plt.ylabel(r'$[\mathrm{HbO}]$')
    plt.xlim([0.05, 0.25])

    # Histograms
    plt.sca(axarr[0,1])
    # plt.hist(transitTimes, bins=20, normed=False)
    plt.hist(pooledTransitTimes, bins=20, normed=False, stacked=True, color=traveledEdgeColors)
    plt.xlabel(r'$\mathrm{Transit\; time}\; [s]$')

    plt.sca(axarr[1,1])
    # plt.hist(finalHb, bins=20, normed=False)
    plt.hist(pooledFinalHb, bins=20, normed=False, stacked=True, color=traveledEdgeColors)
    plt.xlabel(r'$[\mathrm{HbO}]$')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sampledRBC', '-s', help='Name of the sampledRBC directory')
    parser.add_argument('--min-time', type=float, help='Smallest initial RBC time that should be used', \
                        default = 0.0)
    parser.add_argument('--max-time', type=float, help='Largest final RBC time that should be used', \
                        default = 1e100)
    figOptions = FigureOptions(parser)

    args          = parser.parse_args()
    sampledRBCDir = args.sampledRBC
    min_time      = args.min_time
    max_time      = args.max_time
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    sampledRBCsDict = loadSampledRBCs('.', sampledRBCDir)
    displayMoments(sampledRBCsDict, min_time, max_time)
    plotSampledRBCs(sampledRBCsDict, min_time, max_time)

    figOptions.adjustAxes()

    plotName = 'plotTransitTimes_histograms_%s' % sampledRBCDir
    figOptions.saveFig(plotName)
