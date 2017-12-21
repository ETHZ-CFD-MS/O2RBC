#!/usr/bin/env python
"""
Plot longitudinal profiles of PO2. 
"""

import argparse
import matplotlib.pyplot as plt

from HbO2.setup.case import SimulationParametersFactory
from parse.case import time_dirs
from plot.figureoptions import FigureOptions
import plotPO2XProfiles
import plotPO2YProfiles
import plotPO2ProfileAnalytical


def plotPO2SimulatedAndAnalyticalXProfile(caseName, time, radius, **kwargs):
    simParams = SimulationParametersFactory().make_sim_params(caseName)
    fieldName = kwargs.get('fieldName', 'PO2Mean')
    plotPO2XProfiles.plotPO2ProfileEulerian(caseName, float(time), fieldName)
    plotPO2ProfileAnalytical.plotPO2Profile(simParams, radius, **kwargs)


def plotPO2SimulatedAndAnalyticalYProfile(caseName, time, x, **kwargs):
    simParams = SimulationParametersFactory().make_sim_params(caseName)
    fieldNames = kwargs.get('fieldNames', ['PO2'])
    plotPO2YProfiles.plotRadialProfile(caseName, x, float(time), fieldNames, plotFields=['PO2Mean'])
    # plotPO2YProfiles.plotRadialProfileFromPath('PO2RBCAxialAveraged.txt')
    plotPO2ProfileAnalytical.plotPO2YProfile(simParams, x, **kwargs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fields', '-f', nargs='+', help='Fields to plot', default=['PO2'])
    parser.add_argument('--longitudinal', help='Whether to plot a longitudinal profile',
                        action='store_true')
    parser.add_argument('-x', type=float, help='Axial position (for radial profiles)',
                        default=50e-6)
    parser.add_argument('-r', type=float, help='Radial position (for longitudinal profiles)',
                        default=17.6e-6)
    parser.add_argument('--time', '-t', type=float, help='Time for plotting')
    parser.add_argument('--alltimes', '-a', help='Whether to plot all time folders (precedes -t)',
                        action='store_true')
    figOptions = FigureOptions(parser)

    args = parser.parse_args()

    radial_profile = not args.longitudinal
    fieldNames = args.fields
    x = args.x
    r = args.r
    plot_time = args.time
    all_times = args.alltimes
    figOptions.parseOptions(args)

    analyticalLineStyle = {'color': 'k',
                           'linestyle': '-',
                           'dashes': (6, 5)}

    if all_times:
        times = [float(t) for t in time_dirs('domain')]
        if 0 in times:
            times.remove(0)
    else:
        times = [plot_time]

    for time in times:
        if radial_profile:
            plotPO2YProfiles.createPlotDirectory()
            plotPO2SimulatedAndAnalyticalYProfile('.', time, x,
                                                  style=analyticalLineStyle, fieldNames=fieldNames)
            plotName = plotPO2YProfiles.plotName('.', time, x)
        else:
            plotPO2XProfiles.createPlotDirectory()
            plotPO2SimulatedAndAnalyticalXProfile('.', time, r,
                                                  style=analyticalLineStyle, fieldNames=fieldNames)
            plotName = plotPO2XProfiles.plotName('.', time)

        figOptions.adjustAxes()
        figOptions.saveFig(plotName)
        plt.clf()
