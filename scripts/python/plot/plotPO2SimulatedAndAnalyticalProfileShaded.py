#!/usr/bin/env python
"""
Plot profiles of PO2. This script produces only one plot, but with multiple
numerical profiles on top of each other.
"""

import argparse

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import plotPO2ProfileAnalytical
import plotPO2YProfiles
from HbO2.parse import LDRBCVelocityOptions
from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric
from plot.figureoptions import FigureOptions

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fields', '-f', nargs='+', help='Fields to plot', default=['PO2'])
    parser.add_argument('-x', type=float, help='Axial position (for radial profiles)', default=50e-6)
    parser.add_argument('--minTime', type=float, help='Smallest time that should be used',
                        default = 0.0)
    parser.add_argument('--maxTime', type=float, help='Largest time that should be used',
                        default = 2.0)
    parser.add_argument('--step', type=float, help='Time difference between plots that contains all lines',
                        default = 0.1)
    LDRBCVelocityOptions.addOptions(parser)
    figOptions = FigureOptions(parser)

    args = parser.parse_args()
    fieldNames = args.fields
    x          = args.x
    minTime    = args.minTime
    maxTime    = args.maxTime
    step       = args.step
    analyticalOptions = LDRBCVelocityOptions.parseOptions(args)
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    plotPO2YProfiles.createPlotDirectory()

    analyticalLineStyle = {'color': 'k',
                           'linestyle': '-',
                           'dashes': (6, 5)}
    if figOptions.blackWhite:
        PO2RBCStyle = {'color': 'k',
                       'linestyle': '-'}
        PO2MeanStyle = {'color': 'k',
                        'linestyle': '-'}
        cdict = {'red':   ((0.0, 1.0, 1.0),
                           (1.0, 0.6, 0.6)),
                 'green': ((0.0, 1.0, 1.0),
                           (1.0, 0.6, 0.6)),
                 'blue':  ((0.0, 1.0, 1.0),
                           (1.0, 0.6, 0.6))}
        shadingCMapName = 'LightGray'
    else:
        PO2RBCStyle = {'color': 'r',
                       'linestyle': '-'}
        PO2MeanStyle = {'color': 'b',
                       'linestyle': '-'}
        cdict = {'red':   ((0.0, 1.0, 1.0),
                           (1.0, 1.0, 1.0)),
                 'green': ((0.0, 1.0, 1.0),
                           (1.0, 0.6, 0.6)),
                 'blue':  ((0.0, 1.0, 1.0),
                           (1.0, 0.6, 0.6))}
        shadingCMapName = 'WhiteRed'

    shadingCMap = LinearSegmentedColormap(shadingCMapName, cdict)
    plt.register_cmap(cmap=shadingCMap)

    simParams = IOHbO2ParametersAxisymmetric('.')
    plotPO2YProfiles.plotShadedRegionBetweenProfiles('domain', x,
                                minTime, maxTime, step, fieldNames,
                                cmap=shadingCMapName)
    plotPO2YProfiles.plotRadialProfile('domain', x, float(maxTime), fieldNames,
                                       plotFields=['PO2Mean'],
                                       PO2MeanStyle=PO2MeanStyle)
    plotPO2YProfiles.plotRadialProfileFromPath('PO2RBCAxialAveraged.txt',
                                               style=PO2RBCStyle)
    plotPO2ProfileAnalytical.plotPO2YProfile(simParams, x, 
                                             style=analyticalLineStyle,
                                             K0=4.98e6,
                                             convO2Transport=False)

    plt.ylim([35, 70])
    plotName = plotPO2YProfiles.plotName('domain', minTime, x) + '_shaded'

    figOptions.adjustAxes()
    figOptions.saveFig(plotName)
