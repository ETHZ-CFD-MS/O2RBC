#!/usr/bin/env python
"""
Plot data produced by an OpenFOAM probe.
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.interpolate import interp1d

from HbO2.plot.labels import setXLabel, setYLabel
from plot.figureoptions import FigureOptions
from postprocessing.probes import ProbeValue
from postprocessing.probeUtils import getRBCPassingTimes
from utilities.arguments import add_case_argument, add_min_max_time_argument


def plot_probe(probe, probe_idx, min_time, max_time, **kwargs):
    linestyle = kwargs.get('linestyle', {'linestyle': '-'})
    times = probe.times(min_time=min_time, max_time=max_time)
    ymax = 0
    if probe_idx is None:
        probe_idx = probe.position_indices()
    for i in probe_idx:
        values = probe.values(position=i, min_time=min_time, max_time=max_time)
        plt.plot(times, values, **linestyle)
        setXLabel('t', 's')
        setYLabel('PO2', 'mmHg')
        ymax = np.max([ymax, np.max(values)])
    plt.xlim([min(times) - 1e-6, max(times) + 1e-6])
    # plt.ylim([0, 10+np.around(ymax, decimals=-1)])


def plot_passing_RBCs(probe, probe_idx, min_time, max_time, **kwargs):
    RBCPositionsPath = os.path.join(probe.case, 'domain', 'RBCPositions.txt')
    times = probe.times(min_time=min_time, max_time=max_time)
    values = probe.values(position=probe_idx, min_time=min_time, max_time=max_time)
    f = interp1d(times, values, kind='linear')
    front_times = getRBCPassingTimes(probe.positions(probe_idx),
                                     RBCPositionsPath,
                                     'front')
    rear_times  = getRBCPassingTimes(probe.positions(probe_idx),
                                     RBCPositionsPath,
                                     'rear')
    hollow_times  = getRBCPassingTimes(probe.positions(probe_idx),
                                       RBCPositionsPath,
                                       'hollow')
    plt.plot(front_times, f(front_times), 'k.', markersize=10)
    plt.plot(rear_times, f(rear_times),  'k.', markersize=10)
    plt.plot(hollow_times,f(hollow_times),  'r.', markersize=10)

    plot_EAT_bar = False
    if plot_EAT_bar:
        plt.bar(0.038, 29, width=0.002, bottom=26.7, color='0.6')
        plt.text(0.030, 65, 'measured EAT', rotation='horizontal',
                 size=15)
        plt.gca().annotate("",
                           xy=(0.037, 55), xycoords='data',
                           xytext=(0.036, 64), textcoords='data',
                           arrowprops=dict(arrowstyle="fancy",
                                           color="0.6",
                                           shrinkB=5,
                                           connectionstyle="arc3,rad=0.4",
                                           ),
                           )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_case_argument(parser)
    parser.add_argument('--field', '-f', help='Field to plot', default='PO2')
    parser.add_argument('--probeName', '-p', help='Name of the probe directory')
    parser.add_argument('--positions', type=int, nargs='+',
                        help='Positions to plot', default=None)
    add_min_max_time_argument(parser)
    figOptions = FigureOptions(parser)

    args = parser.parse_args()
    case_name = args.caseName
    field_name = args.field
    probe_name = args.probeName
    probe_ids = args.positions
    min_time = args.minTime
    max_time = args.maxTime
    figOptions.parseOptions(args)

    probe = ProbeValue(case_name, probe_name, field_name)
    plot_probe(probe, probe_ids, min_time, max_time)
    # plot_passing_RBCs(probe, probe_ids, min_time, max_time)

    figOptions.adjustAxes()

    plotName = 'plot_{:s}_{:s}'.format(probe_name, field_name)
    figOptions.saveFig(plotName)

