#!/usr/bin/env python
#
# Plot the EATs at a probe for an Eulerian simulation.
#

import argparse
import matplotlib.pyplot as plt
import numpy as np

from plot.figureoptions import FigureOptions
from postprocessing.extractProbeEATs import extractEATs
from postprocessing.loadProbesEuler import loadProbes


def plotProbeEATs(EAT_dict, from_time):

    from_idx = np.searchsorted(EAT_dict['t_front'], from_time)

    plt.scatter(EAT_dict['LD'][from_idx:], EAT_dict['EAT'][from_idx:])
    plt.xlabel('line density [-]')
    plt.ylabel('EAT amplitude [mmHg]')

    # plt.xlim([0, 0.8])
    plt.xlim([0.1, 0.6001])
    plt.ylim([0, 40])

    plt.grid(True)


def plotProbeFlowEATs(EAT_dict, U_RBC, L_RBC, from_time):

    from_idx = np.searchsorted(EAT_dict['t_front'], from_time)

    LD = EAT_dict['LD'][from_idx:]
    flow = [l*U_RBC/L_RBC for l in LD]

    plt.scatter(flow, EAT_dict['EAT'][from_idx:])
    plt.xlabel('RBC flow [cells/s]')
    plt.ylabel('EAT amplitude [mmHg]')

    # plt.xlim([0, 0.8])
    plt.ylim([0, 65])

    plt.grid(True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--field', '-f', help='Field to plot', default='PO2')
    parser.add_argument('--probeName', '-p', help='Name of the probe directory', 
                        default='probeMidstreamPO2')
    parser.add_argument('--from_time', type=float, help='Time from which to plot EATs',
                    default=0.0)
    figOptions = FigureOptions(parser)

    args            = parser.parse_args()
    fieldName       = args.field
    probeName       = args.probeName
    from_time       = args.from_time
    figOptions.parseOptions(args)

    probes = loadProbes('.', probeName, fieldName)
    EAT_dict = extractEATs('.', probes, 0)

    plotProbeEATs(EAT_dict, from_time)
    figOptions.adjustAxes()

    plotName = '%sEATs' % probeName 
    figOptions.saveFig(plotName)

