#!/usr/bin/env python
#
# Plot next to another: longitudinal gradients of different
# simulations
#

import argparse
import matplotlib.pyplot as plt

from plot.figureoptions import FigureOptions
from plotAveragedProbes import plotAveragedProbes, plotAveragedProbePulses
from postprocessing.extractProbeEATs import extractEATs
from postprocessing.extractProbeMinMax import extractProbeMinMax
from postprocessing.loadProbesEuler import loadProbes
from postprocessing import probeUtils


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--field', '-f', help='Field to plot', default='PO2')
    parser.add_argument('--probeNames', '-p', help='Name of the probes to compare', 
                        nargs='+')
    parser.add_argument('--allProbes', '-a', help='Plot all probes (precedes --probeNames)', 
                        action='store_true')
    parser.add_argument('--from_time', type=float, help='Time from which to plot results',
                    default=0.0)
    parser.add_argument('--nearWallPulses', help='Plot near wall pulses', action='store_true')
    parser.add_argument('--cases', help='Name of the cases to compare', 
                        nargs='+')
    figOptions = FigureOptions(parser)

    args      = parser.parse_args()
    fieldName = args.field
    probeNames= args.probeNames
    allProbes = args.allProbes
    from_time = args.from_time
    nearWallPulses  = args.nearWallPulses
    cases     = args.cases
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    linestyle = {'linewidth': 1,
                 'markerfacecolor': 'None',
                 'markeredgecolor': 'None'}
    linecolors = ['k', 'b']
    n = len(cases)
    f, axs = plt.subplots(1,n, sharey=True)

    for i, (case, ax, color) in enumerate(zip(cases, axs, linecolors)):
        plt.sca(ax)
        domainPath = '%s/domain' % case
        
        if allProbes:
            probeNames = probeUtils.probeNames(domainPath, suffix='PO2')

        probes = [loadProbes(domainPath, probeName, fieldName) for probeName in probeNames]
        EAT_dicts = [extractEATs(domainPath, probe, 0) for probe in probes]

        # Plot PO2 max, PO2 mean, inter-RBC PO2 and EATs
        linestyle['color'] = color
        plotAveragedProbes(EAT_dicts, probes, from_time, linestyle=linestyle)

        # Plot amplitude of near-wall pulses
        if nearWallPulses:
            minMaxDicts = [extractProbeMinMax('.', probe, probeIdx=2,
                                              from_time=from_time) for probe in probes]
            plotAveragedProbePulses(minMaxDicts, probes)
        
        if i > 0:
            ax[i].set_ylabel('')

    figOptions.adjustAxes()

    # put A and B at the top left corner of each plot
    axs[0].annotate(r'$\mathbf{A}$', xy=(0.04, 0.91), xycoords='axes fraction')
    axs[1].annotate(r'$\mathbf{B}$', xy=(0.05, 0.91), xycoords='axes fraction')

    plotName = 'compareLongGradient'
    stringList = cases
    stringList.insert(0, plotName)
    plotName = '-'.join(stringList)
    figOptions.saveFig(plotName)





