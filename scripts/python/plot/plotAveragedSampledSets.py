#!/usr/bin/env python
#
# Plots an averaged sampled set.

import argparse
import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from plot.figureoptions import FigureOptions
from postprocessing.extractSampledSetStats import timeAveragedSampledSets
from postprocessing.loadSampledSetFiles import loadSampledSets


## Plot the time-averaged value of a sampledSet
#  @param sampledSetsDict Dictionary with sampled set values
#  @param min_time Time from which data should be averaged
#  @param max_time Time until which data should be averaged
#  @param linestyle Line style
def plotAveragedSampledSet(sampledSetsDict, min_time, max_time, linestyle={'linestyle':'-'}):

    positions = sampledSetsDict['positions']
    plotPositions = np.linspace(positions[0], positions[-1], 300)
    averagedValues = timeAveragedSampledSets(sampledSetsDict, min_time, max_time)
    # print averagedValues

    f = interp1d(positions[:], averagedValues, kind='cubic')
    plt.plot(1e6*plotPositions, f(plotPositions), **linestyle)

    plt.xlabel(r'$y \; [\mu \mathrm{m}]$')
    plt.ylabel(r'$\mathrm{PO}_2 \; [\mathrm{mm\,Hg}]$')

    # plt.ylim([0, max(averagedValues)])
    # plt.ylim([-10+np.around(min(averagedValues), decimals=-1), 10+np.around(max(averagedValues), decimals=-1)])
    plt.ylim([20, 60])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sampledSet', '-s', help='Name of the sampledSet directory')
    parser.add_argument('--fileName', '-f', help='Name of the sampledSet files')
    parser.add_argument('--min-time', type=float, help='Smallest time that should be used',
                        default = 0.0)
    parser.add_argument('--max-time', type=float, help='Largest time that should be used',
                        default = 1e100)
    parser.add_argument('--cases', help='Name of the cases to plot',
                        nargs='+', default=['.'])
    figOptions = FigureOptions(parser)

    args          = parser.parse_args()
    sampledSetDir = args.sampledSet
    fileName      = args.fileName
    min_time      = args.min_time
    max_time      = args.max_time
    cases         = args.cases
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    styles = []
    styles.append({'color': 'b',
                  'linestyle': '-',})
    styles.append({'color': 'r',
                  'linestyle': '-',
                  'dashes': (4,2.5)})
    styles.append({'color': 'k',
                  'linestyle': '-',
                  'dashes': (1,2)})
    styles.append({'color': 'g',
                  'linestyle': '-.',
                  'dashes': (4,2,1,2)})

    for case, style in zip(cases, styles):
        domainPath = '%s/domain' % case
        sampledSetsDict = loadSampledSets(domainPath, sampledSetDir, fileName)
        plotAveragedSampledSet(sampledSetsDict, min_time, max_time, style)

    # plt.ylim([20, 45])
    figOptions.adjustAxes()
    figOptions.setGrid()

    plotName = 'plotAveraged_%s' % fileName
    figOptions.saveFig(plotName)


