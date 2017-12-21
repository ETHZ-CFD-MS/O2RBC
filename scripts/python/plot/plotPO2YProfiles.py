#!/usr/bin/env python
"""
Plot radial profiles of PO2. 
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os

from HbO2.plot import labels
from parse.case import time_dirs
from parse.readfile import load_numeric_to_dictionary
from plot.figureoptions import FigureOptions
from postprocessing import readSampleFiles

plotDirName = 'pySamplePlots_y_axis'
setDirectory = 'sets_y_axis'


def plotRadialProfile(caseName, x, time, fieldNames, plotFields='all', **kwargs):
    if plotFields == 'all':
        plotFields = fieldNames

    # extract sample data
    xString = '%i' % (1e6*x)
    setNames = ['x_%s' % xString]
    labelText = [r'x = %s {\large $\mu$}m' % xString]
    PO2Style = kwargs.get('PO2Style', {'color': (1, 0.5, 0.5),
                                       'linestyle': '-',
                                       'linewidth': 0.3,
                                       'alpha': 0.1})
    PO2MeanStyle = kwargs.get('PO2MeanStyle', {'color': 'b',
                                               'linestyle': '-'})

    for setName, label in zip(setNames, labelText):
        sample_data = readSampleFiles.loadSampleFile(caseName, setDirectory,
                            time, setName, fieldNames)

        y_coords = sample_data[:,0]
        data     = sample_data[:,1:]

        # convert to microns
        y_coords = y_coords*1e6
        # plot extracted profiles
        for i, fieldName in enumerate(fieldNames):
            if fieldName == 'PO2' and 'PO2' in plotFields:
                plt.plot(y_coords, data[:,i], label=label, **PO2Style)
            elif fieldName == 'PO2Mean' and 'PO2Mean' in plotFields:
                yMin = 2.6
                # yMin = 0.0
                filteredCoords = [y          for y in y_coords if y >= yMin]
                filteredValues = [data[j, i] for j, y in enumerate(y_coords) \
                                                      if y >= yMin]
                plt.plot(filteredCoords, filteredValues, label=label, 
                         **PO2MeanStyle)

    labels.setXLabel('y', 'um')
    labels.setYLabel('PO2', 'mmHg')


def plotRadialProfileFromPath(pathToData, **kwargs):
    style = kwargs.get('style', {'color': 'r',
                                 'linestyle': '-'})
    data_dict = load_numeric_to_dictionary(open(pathToData))
    x = 1e6*data_dict['r']
    vals = data_dict['PO2']
    plt.plot(x, vals, **style)


def plotShadedRegionBetweenProfiles(caseName, x, timeStart, timeEnd, \
                                     step, fieldNames, cmap='Reds'):
    nPO2Bins = 500
    nLevels  = 50
    nTimeStep = int(np.round(1. + (timeEnd - timeStart)/step))
    times = np.linspace(timeStart, timeEnd, nTimeStep)
    dataDict = readSampleFiles.loadSampleFileFromTimes(caseName, setDirectory, \
                                                times, setName(x), fieldNames)
    minPO2  = np.min(dataDict['PO2'])
    maxPO2  = np.max(dataDict['PO2'])
    PO2BinEdges   = np.linspace(minPO2, maxPO2, nPO2Bins+1)
    PO2BinCenters = 0.5*(PO2BinEdges[0:-1] + PO2BinEdges[1:])
    frequency = countFrequency(dataDict, nPO2Bins)
    frequency = smoothFrequency(frequency)
    shade = np.transpose(mapFrequencyToShade(frequency))
    levels = np.linspace(0, 1, nLevels)
    plt.contourf(1e6*dataDict['x'], PO2BinCenters, shade, levels=levels, 
                 cmap=cmap)

## Return the frequency of the PO2 values in an ny x nPO2Bins numpy array.
#  The frequency is interpolated as follows: if two consecutive samples
#  fall into different bins, all the traversed bins are incremented by
#  1/(# traversed bins)
def countFrequency(dataDict, nPO2Bins):
    ny      = dataDict['PO2'].shape[0]
    nSample = dataDict['PO2'].shape[1]
    minPO2  = np.min(dataDict['PO2'])
    maxPO2  = np.max(dataDict['PO2'])
    freq = np.zeros((ny, nPO2Bins)) 
    binNumber = transformToBinIndex(dataDict['PO2'], minPO2, maxPO2, nPO2Bins)
    for i in range(ny):
        for j in range(nSample-1):
            binRange = absoluteInclusiveRange(binNumber[i, j], binNumber[i, j+1])
            for binI in binRange:
                freq[i, int(binI)] += 1./len(binRange)
    return freq

## Smooth the frequency on vertical slices. Smoothing is only applied to 
#  the nonzero portion of the signal and when its length is larger than
#  the smoothing window length.
def smoothFrequency(freq):
    windowLen = 21
    w = np.hanning(windowLen)
    for i in range(freq.shape[0]):
        nzRange = nonZeroIndexRange(freq[i,:])
        if len(nzRange) > windowLen:
            x = freq[i,nzRange]
            s = np.r_[x[windowLen-1:0:-1],x,x[-1:-windowLen:-1]]
            smoothedSlice = np.convolve(w/w.sum(), s, mode='valid')
            freq[i,nzRange] = smoothedSlice[windowLen/2-1:-(windowLen/2+1)]
    return freq

# Return the integer bin index for the values in array and the bins defined
# by the remaining arguments.
def transformToBinIndex(array, minValue, maxValue, nBins):
    binIndex = np.floor(nBins*(array - minValue)/(maxValue - minValue))
    binIndex[binIndex == nBins] = nBins - 1
    return binIndex.astype(int)

# Return a range between integers that includes both of them.
# The arguments need not be in increasing order.
def absoluteInclusiveRange(a, b):
    if a == b:
        return [a]
    elif a > b:
        return range(b, a+1)
    else:
        return range(a, b)

def mapFrequencyToShade(frequency):
    return frequency/(10 + frequency)

## Return the largest index range so that the first and last elements
#  are nonzero.
def nonZeroIndexRange(x):
    nonZeroIndices = np.nonzero(x[:])[0]
    return np.arange(nonZeroIndices[0],nonZeroIndices[-1]+1)


def setName(x):
    return 'x_%i' % (1e6*x)


def createPlotDirectory():
    plotDirPath = os.path.join('domain', plotDirName)
    if not os.path.exists(plotDirPath):
        os.makedirs(plotDirPath)


def plotName(caseName, time, x):
    return os.path.join(caseName, 'domain', plotDirName,
                        'PO2YProfiles_x_%g_t_%.5f' % (1e6*x, time))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fields', '-f', nargs='+', help='Fields to plot', default=['PO2'])
    parser.add_argument('-x', type=float, help='x-coordinate to plot', default=50e-6)
    parser.add_argument('--alltimes', '-a', help='Whether to plot all time folders',
                        action='store_true')
    parser.add_argument('--min-time', type=float, help='Smallest time that should be used',
                        default = 0.0)
    parser.add_argument('--max-time', type=float, help='Largest time that should be used',
                        default = 2.0)
    parser.add_argument('--step', type=float, help='Time difference between plots that contains all lines', \
                        default = 0.1)
    figOptions = FigureOptions(parser)
    args = parser.parse_args()
    fieldNames = args.fields
    x = args.x
    all_times = args.alltimes
    min_time  = args.min_time
    max_time  = args.max_time
    step = args.step
    figOptions.parseOptions(args)

    createPlotDirectory()

    if all_times:
        times = time_dirs('domain')
        if '0' in times:
            times.remove('0')
    else:
        n = int(round((max_time-min_time)/step) + 1)
        times = np.linspace(min_time, max_time, n)
    for time in times:
        plotRadialProfile('domain', x, float(time), fieldNames)
        figOptions.adjustAxes()
        figOptions.setGrid()
        figOptions.saveFig(plotName('domain', float(time), x))
        plt.clf()


