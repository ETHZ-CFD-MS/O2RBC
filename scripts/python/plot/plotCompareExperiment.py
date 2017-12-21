#!/usr/bin/env python
#
# Compare results for a csv file with simulated results.
# Currently working for PO2 min/max
#

import argparse

import matplotlib.pyplot as plt
import numpy as np

from parse.readPO2CSV import readPO2CSV
from plot.figureoptions import FigureOptions
from plot.plotProbeMinMax import plotProbeMinMax
from postprocessing.extractProbeEATs import extractEATs
from postprocessing.loadProbesEuler import loadProbes


## Bar plot of PO2 and inter-RBC PO2 averaged of hematocrit bins
#  @param csv_dict Dictionary typically produced by readPO2CSV. Entries are 1D numpy arrays.
def plotAveragedBar(csv_dict):

    Ht = csv_dict['Hematocrit']
    RBC_PO2 = csv_dict['RBC']
    interRBC_PO2 = csv_dict['interRBC']
    Ht_bins = [5, 20, 40, 63]

    (RBC_PO2_av, RBC_PO2_std)           = binnedAverage(Ht, RBC_PO2, Ht_bins)
    (interRBC_PO2_av, interRBC_PO2_std) = binnedAverage(Ht, interRBC_PO2, Ht_bins)
    EAT = [RBC_PO2_av[i] - interRBC_PO2_av[i] for i in range(len(RBC_PO2_av))]

    ax = plt.gca()
    ind = np.arange(len(Ht_bins)-1)
    width = 0.6
    p1 = ax.bar(ind, interRBC_PO2_av, width, color=(0.0,0.0,1.0,0.6), label='inter-RBC PO2')
    p2 = ax.bar(ind, EAT, width, color='1', bottom=interRBC_PO2_av, yerr=RBC_PO2_std, label='RBC PO2')
    ax.errorbar(ind+0.5*width, RBC_PO2_av, yerr=RBC_PO2_std, fmt=None, ecolor='k')

    ax.set_xticks(ind+0.5*width)
    ax.set_xticklabels( (r'$0.05-0.2$', r'$0.2-0.4$', r'$0.4-0.6$') )
    ax.set_xlim(min(ind)-0.3*width, max(ind)+1.3*width)
    plt.xlabel(r'$\mathrm{linear\; density\; [-]}$')
    ax.tick_params(\
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='on') # labels along the bottom edge are off


    # ax.legend(['a', 'b'], loc='upper left')
    # ax.legend(loc='upper left')

## Average ydata in a numpy array along bins for corresponding xdata
#  @param xdata x-data (numpy array)
#  @param ydata y-data (numpy array)
#  @param bins Bins for x-data
#  @return (averages,stdevs) Tuple with averages and std. deviations of y-data (lists)
def binnedAverage(xdata, ydata, bins):
    averages = []
    stdevs = []
    for i in range(len(bins)-1):
        idx_bin = np.where(np.logical_and(xdata >= bins[i], xdata < bins[i+1]))
        averages.append(np.mean(ydata[idx_bin]))
        stdevs.append(np.std(ydata[idx_bin]))

    return (averages, stdevs)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--field', '-f', help='Field to plot', default='PO2')
    parser.add_argument('--probeName', '-p', help='Name of the probe to plot')
    parser.add_argument('--from_time', type=float, help='Time from which to plot results',
                    default=0.0)
    figOptions = FigureOptions(parser)

    args            = parser.parse_args()
    fieldName       = args.field
    probeName       = args.probeName
    from_time       = args.from_time
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    probes = loadProbes('domain', probeName, fieldName)
    EAT_dict = extractEATs('domain', probes, 0)
    f, (ax1, ax2) = plt.subplots(1,2, sharey=True)
    plt.sca(ax1) # set current axes
    plotType = 'minMaxTissue'
    # plotType = 'minMax'
    plotProbeMinMax(EAT_dict, None, from_time, plotType)
    ax1.set_xlim([0.15, 0.6])

    plt.sca(ax2) # set current axes
    fileName = 'flow_hematocrit_EAT.csv'
    csv_dict = readPO2CSV(fileName, delimiter=',')
    plotAveragedBar(csv_dict)

    figOptions.adjustAxes()
    figOptions.setGrid()

    # put (a) and (b) at the top left corner of each plot
    ax1.annotate(r'$\mathbf{A}$', xy=(0.03, 0.91), xycoords='axes fraction')
    ax2.annotate(r'$\mathbf{B}$', xy=(0.03, 0.91), xycoords='axes fraction')

    plotName = 'compareExperimentWith%sMinMax' % probeName
    figOptions.saveFig(plotName)



