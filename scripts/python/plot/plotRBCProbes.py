#!/usr/bin/env python
#
# Plot the values from the function object sampleRBCField as a function of time
#

import argparse

import matplotlib.pyplot as plt

from extractRBCProbeStats import RBCProbeValues
from loadRBCProbes import loadRBCProbes
from plot.figureoptions import FigureOptions


def plotRBCProbe(RBCProbesDict, edgeIDs, sCoord, min_time, max_time):
    legends = []
    for eI in edgeIDs:
        probeDict = RBCProbeValues(RBCProbesDict, eI, sCoord)
        times      = probeDict['times']
        meanValues = probeDict['meanValues']

        plotTimeIdx = [i for i,t in enumerate(times) if t >= min_time and t <= max_time]
        if len(plotTimeIdx) == 0:
            raise ValueError('No time selected, check the arguments min_time and max_time.')
        timesPlot = times[plotTimeIdx]
        plt.plot(timesPlot, meanValues[plotTimeIdx], '-o')
        legends.append('edge %i' % eI)

    plt.xlim([min_time - 1e-6, max_time + 1e-6])
    plt.ylim([0.2, 0.7])
    plt.legend(legends, loc='lower right')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--field', '-f', help='Field to plot', default='Hb')
    parser.add_argument('--probeName', '-p', help='Name of the probe directory', 
                        default='probe05PO2')
    parser.add_argument('--edgeIndex', '-e', type=int, nargs='+',
                         help='Edge indices of the probes to plot')
    parser.add_argument('--sCoord', '-s', nargs='+',
                        help='Curvilinear coordinates of the probes to sample, or ''any''',
                        default='any')
    parser.add_argument('--min-time', type=float, help='Smallest time that should be used', \
                        default = 0.0)
    parser.add_argument('--max-time', type=float, help='Largest time that should be used', \
                        default = 1.0)
    figOptions = FigureOptions(parser)

    args = parser.parse_args()
    fieldName = args.field
    probeName = args.probeName
    edgeIDs = args.edgeIndex
    try:
        sCoord = float(args.sCoord)
    except ValueError:
        sCoord = args.sCoord
    min_time = args.min_time
    max_time = args.max_time

    figOptions.parseOptions(args)
    figOptions.applyOptions()

    RBCProbesDict = loadRBCProbes('.', probeName, fieldName)
    plotRBCProbe(RBCProbesDict, edgeIDs, sCoord, min_time, max_time)

    figOptions.adjustAxes()
    figOptions.setGrid()
    plotName = '%sPlot' % probeName
    figOptions.saveFig(plotName)

