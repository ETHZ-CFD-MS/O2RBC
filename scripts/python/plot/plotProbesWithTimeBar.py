#!/usr/bin/env python
"""
Plot data produced by an OpenFOAM probe for all time steps, with a
vertical line that indicates the current time.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np

from postprocessing.probes import ProbeValue
from plot.figureoptions import FigureOptions
from plot.plotProbe import plot_probe
from utilities.arguments import add_case_argument, add_min_max_time_argument


parser = argparse.ArgumentParser()
add_case_argument(parser)
parser.add_argument('--field', '-f', help='Field to plot', default='PO2')
parser.add_argument('--probeName', '-p', help='Name of the probe directory', 
                    default='probeMidstreamPO2')
parser.add_argument('--positions', type=int, nargs='+', 
                    help='Positions to plot', default=None)
add_min_max_time_argument(parser)
parser.add_argument('--step', type=float, help='Time difference between plots that contains all lines',
                default = 0.01)
figOptions = FigureOptions(parser)

args = parser.parse_args()
case_name = args.caseName
field_name = args.field
probe_name = args.probeName
probe_ids = args.positions
min_time = args.minTime
max_time = args.maxTime
step = args.step

figOptions.parseOptions(args)

n = int(round((max_time-min_time)/step) + 1)
times = np.linspace(min_time, max_time, n)
for i in xrange(n):
    t = times[i]
    print t
    plt.cla()
    probe = ProbeValue(case_name, probe_name, field_name)
    plot_probe(probe, probe_ids, min_time, max_time)
    plt.plot([t, t], [0 , 70], color='0.4')
    if i == 0:
        figOptions.adjustAxes()

    plotName = '%sPlot' % probe_name
    figOptions.saveFig(plotName)

    newPlotName = '%s_%04i.png' % (probe_name, i)
    os.system('mv %s.png %s' % (plotName, newPlotName))
    print 'Moved plot to %s' % newPlotName

