#!/usr/bin/env python
"""
Plot data produced by an OpenFOAM probe.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt

from HbO2.plot import labels
from loadPostProcessingFiles import loadSampleFile
from plot.figureoptions import FigureOptions


def plotProfile(caseName, time, profileName, linestyle={'linestyle':'-'}):
    ## user-editable part:
    # Indicate which coordinate should be used for plotting 
    # If the three coordinates are written: 0 -> x, 1 -> y, 2 -> z
    # If only one coordinate was written, use 0
    coordIdx = 1

    labels.setXLabel('x', 'um')
    if coordIdx == 1:
        labels.setXLabel('y', 'um')
    labels.setYLabel('PO2', 'mmHg')

    profileData = loadSampleFile(caseName, time, profileName)
    coords = profileData[:,0:-1]
    values = profileData[:,-1]

    plt.plot(1e6*coords[:,coordIdx], values, **linestyle)
    plt.ylim([-5+np.around(min(values), decimals=-1),
               5+np.around(max(values), decimals=-1)])


if __name__ == '__main__':
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', '-p', help='Name of the profile file')
    parser.add_argument('--time', '-t', help='Time for which to plot the profile')
    figOptions = FigureOptions(parser)

    args        = parser.parse_args()
    profileName = args.profile
    time        = args.time
    figOptions.parseOptions(args)

    figOptions.applyOptions()

    caseName = '.'

    plotProfile(caseName, time, profileName)

    figOptions.adjustAxes()
    plotName = '%sPlot' % profileName
    figOptions.saveFig(plotName)
