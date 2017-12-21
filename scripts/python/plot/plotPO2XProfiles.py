#!/usr/bin/env python
#
# Plot longitudinal profiles of PO2. 
#

import argparse
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
import os
from scipy.interpolate import interp1d

from HbO2.plot import labels
from parse.case import time_dirs
from plot.figureoptions import FigureOptions
from postprocessing import profileUtils, readSampleFiles

plotDirName = 'pySamplePlots_x_axis'


def plotPO2ProfileEulerian(caseName, time, fieldName, **kwargs):
    defaultSetNames = ['centerline', 'nearWall', 'fiveMuFromWall', 'tenMuFromWall']
    defaultLabels   = ['centerline', r'1 $\muup$m from wall',
                       r'5 $\muup$m from wall', r'10 $\muup$m from wall']
    defaultStyles   = []
    defaultStyles.append({'color': 'b',
                          'linestyle': '-',})
    defaultStyles.append({'color': 'r',
                          'linestyle': '-',
                          'dashes': (4,2.5)})
    defaultStyles.append({'color': 'g',
                          'linestyle': '-.',
                          'dashes': (4,2,1,2)})
    defaultStyles.append({'color': 'k',
                          'linestyle': '-',
                          'dashes': (1,2)})
    alpha = kwargs.get("alpha", 1)
    setNames = kwargs.get("setNames", defaultSetNames)
    labelNames   = kwargs.get("labelNames",   defaultLabels)
    styles   = kwargs.get("styles",   defaultStyles)

    # extract sample data
    setDirectory = 'sets_x_axis'
    # setNames     = ['centerline', 'nearWall', 'fiveMuFromWall', 'tenMuFromWall', 'tissueExterior']
    # labelNames       = ['centerline', r'1 {\small $\mu$}m from wall', r'5 {\small $\mu$}m from wall', \
                    # r'10 {\small $\mu$}m from wall', r'15 {\small $\mu$}m from wall']
    # setNames     = ['centerline', 'nearWall', 'fiveMuFromWall', 'tenMuFromWall']
    # labelNames       = ['centerline', r'1 {\large $\mu$}m from wall',
                    # r'5 {\large $\mu$}m from wall', r'10 {\large $\mu$}m from wall']
    # setNames     = ['tenMuFromWall']
    # labelNames       = [r'15 {\large $\mu$}m from wall']

    ymax = 0

    for setName, label, style in zip(setNames, labelNames, styles):
        sample_data = readSampleFiles.loadSampleFile(caseName, setDirectory, time, setName, [fieldName])

        x_coords = sample_data[:,0]
        PO2      = sample_data[:,1]
        ymax     = np.max([ymax, np.max(PO2)])

        # convert to microns
        x_coords = x_coords*1e6

        # plot extracted profiles
        plt.plot(x_coords, PO2, label=label, linewidth=1, alpha=alpha, **style)

        if setName == 'centerline':
            # interpolate profile
            f = interp1d(x_coords, PO2, kind='cubic', bounds_error=False,
                    fill_value=0.0)
            # plotRBCPositions(caseName, time, 'front', f, color='k')
            # plotRBCPositions(caseName, time, 'rear', f, color='k')
            # plotRBCPositions(caseName, time, 'hollow', f, color='r')
            # plotRBCThickLines(caseName, time)

    xmin = x_coords[0]
    xmax = x_coords[-1]

    labels.setXLabel('x', 'um')
    labels.setYLabel('PO2', 'mmHg')
    plt.xlim([xmin, xmax])
    plt.ylim([0, 10+np.around(ymax, decimals=-1)])
    create_legend(plt.gca(), bbox_to_anchor=(0.98, 0.96))


def plotRBCPositions(caseName, time, RBCPart, f, color='k'):
    RBCPositionsPath = '%s/RBCPositions.txt' % caseName
    RBCPositions = profileUtils.getRBCPositions(time,
                RBCPositionsPath, RBCPart) 
    RBCPositions = [1e6*x for x in RBCPositions]
    plt.plot(RBCPositions, f(RBCPositions), '.', color=color, markersize=8)


def plotRBCThickLines(caseName, time):
    # TODO: repair this
    RBCPositionsPath = '%s/RBCPositions.txt' % caseName
    ax = plt.gca()
    RBCFronts = profileUtils.getRBCPositions(time,
                RBCPositionsPath, 'front') 
    RBCRears = profileUtils.getRBCPositions(time,
                RBCPositionsPath, 'rear') 
    RBCFronts = 1e6*RBCFronts
    RBCRears  = 1e6*RBCRears
    for i in range(len(RBCFronts)):
        rect = patches.Rectangle((RBCRears[i], 0), RBCFronts[i]-RBCRears[i], 3, color='black')
        ax.add_patch(rect)


def createPlotDirectory():
    plotDirPath = os.path.join('domain', plotDirName)
    if not os.path.exists(plotDirPath):
        os.makedirs(plotDirPath)


def plotName(caseName, time):
    return os.path.join(caseName, 'domain', plotDirName,
                        'PO2XProfiles_t_{:.3f}'.format(time))


def create_legend(ax, **kwargs):
    """
    Create a legend for the axis in argument.

    Args:
        ax (pyplot.Axis): axis object

    Returns:
        created legend object
    """
    fontP = FontProperties()
    fontP.set_size('x-small')
    props = dict({'prop': fontP,
                  'bbox_to_anchor': (1, 0.80),
                  'fancybox': True,
                  'borderpad': 0.2,
                  'labelspacing': 0.3,
                  'handlelength': 3.6},
                 **kwargs)
    legend = ax.legend(**props)
    legend.get_frame().set_linewidth(0.5)
    return legend


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--field', '-f', help='Field to plot', default='PO2')
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
    fieldName = args.field
    all_times = args.alltimes
    min_time  = args.min_time
    max_time  = args.max_time
    step      = args.step
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
        plotPO2ProfileEulerian('domain', float(time), fieldName)
        figOptions.adjustAxes()
        figOptions.saveFig(plotName('.', float(time)))
        plt.clf()

