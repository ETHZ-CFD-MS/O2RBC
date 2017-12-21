#!/usr/bin/env python
#
# Plot next to another: x-profiles (snapshot) and probe comparison (average)
#

import argparse
import matplotlib.pyplot as plt

from plot.figureoptions import FigureOptions
from plot.plotPO2XProfiles import plotPO2ProfileEulerian
from postprocessing.extractProbeEATs import extractEATs
from postprocessing.loadProbesEuler import loadProbes
from postprocessing import probeUtils


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--time', type=float, help='Time for the snapshot')
    parser.add_argument('--from_time', type=float, help='Time from which to average results', 
                    default=0.0)
    figOptions = FigureOptions(parser)

    args            = parser.parse_args()
    time            = args.time
    from_time       = args.from_time
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    f, (ax1, ax2) = plt.subplots(1,2, sharey=True)
    plt.sca(ax1) # set current axes
    plotPO2ProfileEulerian('.', time)

    plt.sca(ax2) # set current axes
    probeNames = probeUtils.probeNames('.', suffix='PO2')
    fieldName = 'PO2'
    probes = [loadProbes('.', probeName, fieldName) for probeName in probeNames]
    EAT_dicts = [extractEATs('.', probe, 0) for probe in probes]
    plotCompareProbes(EAT_dicts, probes, from_time)
    ax2.set_ylabel('')

    figOptions.adjustAxes()
    figOptions.setGrid()

    # put (a) and (b) at the top left corner of each plot
    ax1.annotate(r'$(a)$', xy=(0.03, 0.91), xycoords='axes fraction')
    ax2.annotate(r'$(b)$', xy=(0.03, 0.91), xycoords='axes fraction')

    plotName = 'PO2XProfilesAndCompareProbes'
    figOptions.saveFig(plotName)




