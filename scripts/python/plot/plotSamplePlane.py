#!/usr/bin/env python
#
# Plot data produced by an OpenFOAM sample plane in raw format
#

import argparse
import matplotlib.pyplot as plt
from matplotlib.delaunay import *
import numpy as np
import scipy

from plot.figureoptions import FigureOptions
from postprocessing.readPostProcessingFile import readPostProcessingFile


# Read a sample plane file produced by OpenFOAM and return its content in a dictionary with
# keys 'header' and 'values'
def loadSamplePlane(caseName, dirName, timeName, planeName):
    filePath = '%s/postProcessing/%s/%s/%s' % (caseName, dirName, timeName, planeName)
    samplePlaneDict = readPostProcessingFile(filePath)
    return samplePlaneDict

def defGrid(x, y, delta):
    xvalues = np.arange(min(x), max(x), delta)
    yvalues = np.arange(min(y), max(y), delta)
    xi, yi = scipy.meshgrid(xvalues, yvalues)
    return xi, yi

def interpolateToMeshGrid(x, y, w, xi, yi):
    tri = Triangulation(x, y)
    interp = tri.nn_interpolator(w)
    wi = interp(xi, yi)
    return wi

def plotSamplePlane(samplePlaneDict):
    x = samplePlaneDict['values'][:,0]
    y = samplePlaneDict['values'][:,1]
    z = samplePlaneDict['values'][:,2]
    w = samplePlaneDict['values'][:,3]

    xi, yi = defGrid(y, z, 1e-6)
    wi     = interpolateToMeshGrid(y, z, w, xi, yi)

    nLevels = 20
    levels = np.arange(min(w), max(w), (max(w) - min(w))/nLevels)
    CS = plt.contour(xi, yi, wi, levels)
    plt.clabel(CS, inline=1, fontsize=10)


if __name__ == '__main__':
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--sampleName', '-s', help='Name of the sample plane directory')
    parser.add_argument('--time', '-t', help='Time directory')
    parser.add_argument('--planeName', '-p', help='File name of the sample plane')
    figOptions = FigureOptions(parser)

    args            = parser.parse_args()
    sampleName      = args.sampleName
    timeName        = args.time
    planeName       = args.planeName
    figOptions.parseOptions(args)

    figOptions.applyOptions()

    caseName = '.'

    samplePlaneDict = loadSamplePlane(caseName, sampleName, timeName, planeName)
    plotSamplePlane(samplePlaneDict)

    figOptions.adjustAxes()
    figOptions.setGrid()

    plotName = '%sPlot' % planeName
    figOptions.saveFig(plotName)


