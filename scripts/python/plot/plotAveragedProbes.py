#!/usr/bin/env python
#
# Plots the comparison between two probes.
# Can be used to assess longitudinal gradients.
#

import argparse

import matplotlib.pyplot as plt
import numpy as np

from HbO2.plot import labels
from postprocessing.extractProbeEATs import extractEATs
from postprocessing.extractProbeMinMax import extractProbeMinMax
from postprocessing.loadProbesEuler import loadProbes
from postprocessing import probeUtils
from plot.figureoptions import FigureOptions


## Plot the values of RBC PO2, inter-RBC PO2 and mean PO2 at multiple probes using bars.
#  @param EAT_dicts List with dictionaries read by extractProbeEATs.
#  @param from_time Time from which data should be plotted.
def plotAveragedProbesBars(EAT_dicts, probeNames, from_time):

    nProbes = len(probeNames)

    from_idx = np.searchsorted(EAT_dicts[0]['t_front'], from_time)

    PO2_max  = [EAT_dict['PO2_max'][from_idx:]  for EAT_dict in EAT_dicts]
    PO2_min  = [EAT_dict['PO2_min'][from_idx:]  for EAT_dict in EAT_dicts]
    PO2_mean = [EAT_dict['PO2_mean'][from_idx:] for EAT_dict in EAT_dicts]

    mean_PO2_max  = [np.mean(P) for P in PO2_max]
    mean_PO2_min  = [np.mean(P) for P in PO2_min]
    mean_PO2_mean = [np.mean(P) for P in PO2_mean]

    std_PO2_max  = [np.std(P) for P in PO2_max]
    std_PO2_min  = [np.std(P) for P in PO2_min]
    std_PO2_mean = [np.std(P) for P in PO2_mean]

    probe_means = [[PO2_max, PO2_min, PO2_mean] \
                    for PO2_max, PO2_min, PO2_mean in zip(mean_PO2_max, mean_PO2_min, mean_PO2_mean)]
    probe_std   = [[PO2_max, PO2_min, PO2_mean] \
                    for PO2_max, PO2_min, PO2_mean in zip(std_PO2_max, std_PO2_min, std_PO2_mean)]
    
    ind = np.arange(3) # three statistics are plotted: max, min, mean
    width = 1./(nProbes+2)

    # http://matplotlib.org/examples/api/barchart_demo.html
    fig, ax = plt.subplots()
    rects = []
    colors = ['%f' % c for c in np.linspace(0,1, nProbes+2)]
    colors = colors[1:]

    for i in range(nProbes):
        rects.append(ax.bar(ind+i*width, probe_means[i], width, color=colors[i], yerr=probe_std[i]))

    # rects1 = ax.bar(ind,       probe_means[0], width, color='k', yerr=probe_std)
    # rects2 = ax.bar(ind+width, probe_means[1], width, color='0.8', yerr=probe_std)

    ax.set_ylabel(r'$\mathrm{PO}_2 \; [\mathrm{mm\,Hg}]$')
    ax.set_xticks(ind+0.5*nProbes*width)
    ax.set_xticklabels( (r'RBC PO$_2$', r'inter-RBC PO$_2$', r'mean PO$_2$') )

    ax.legend( rects, probeNames )

    for rect in rects:
        autolabel(rect)
    
    print 'Mean PO2 max = ', mean_PO2_max
    print 'Std PO2 max  = ', std_PO2_max
    print 'Mean PO2 min = ', mean_PO2_min
    print 'Std PO2 min  = ', std_PO2_min
    print 'Mean PO2 mean = ', mean_PO2_mean
    print 'Std PO2 mean  = ', std_PO2_mean

def autolabel(rects):
    # attach some text labels
    ax = plt.gca()
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%.1f' % height,
                ha='center', va='bottom')

## Plot the values of RBC PO2, inter-RBC PO2 and mean PO2 at multiple probes using curves.
#  @param EAT_dicts List with dictionaries read by extractProbeEATs.
#  @param probes List of dictionaries produces by loadProbesEuler.py
#  @param from_time Time from which data should be plotted.
def plotAveragedProbes(EAT_dicts, probes, from_time, **kwargs):

    nProbes = len(probes)

    # extract probe positions (scaled to mum)
    xProbes = [1e6*probe['positions'][0][0] for probe in probes]

    # extract PO2 max, min and mean
    from_idx = np.searchsorted(EAT_dicts[0]['t_front'], from_time)

    PO2_max  = [EAT_dict['PO2_max'][from_idx:]  for EAT_dict in EAT_dicts]
    PO2_min  = [EAT_dict['PO2_min'][from_idx:]  for EAT_dict in EAT_dicts]
    PO2_mean = [EAT_dict['PO2_mean'][from_idx:] for EAT_dict in EAT_dicts]

    mean_PO2_max  = np.asarray([np.mean(P) for P in PO2_max])
    mean_PO2_min  = np.asarray([np.mean(P) for P in PO2_min])
    mean_PO2_mean = np.asarray([np.mean(P) for P in PO2_mean])

    mean_EAT = mean_PO2_max - mean_PO2_min

    style = kwargs.get('linestyle', {'markerfacecolor': 'None',
                                     'markersize': 5,
                                     'linewidth': 0.4})
    plt.plot(xProbes, mean_PO2_max, 'ro-', label=r'RBC', **style)
    # plt.plot(xProbes, mean_PO2_max, 'r-', label=r'RBC', **style)
    plt.plot(xProbes, mean_PO2_mean,'ks-', label=r'mean', **style)
    plt.plot(xProbes, mean_PO2_min, 'b^-', label=r'inter-RBC', **style)
    # plt.plot(xProbes, mean_PO2_min, 'b-', label=r'inter-RBC', **style)
    plt.plot(xProbes, mean_EAT,     'k--', dashes=(5,3), label='EAT', **style)

    labels.setXLabel('x', 'um')
    labels.setYLabel('PO2', 'mmHg')

    # plt.legend()
    plt.ylim([0, 80])
    # plt.xlim([8.6, 91.4])
    # plt.xlim([10, 90])

    print 'Mean PO2 max = ', mean_PO2_max
    print 'Mean PO2 min = ', mean_PO2_min
    print 'Mean PO2 mean = ', mean_PO2_mean
    print 'Mean EAT      = ', mean_EAT
    print 'Average mean PO2 max = ', np.mean(mean_PO2_max)
    print 'Average mean PO2 min = ', np.mean(mean_PO2_min)
    print 'Average mean PO2 mean = ', np.mean(mean_PO2_mean)
    print 'Average mean EAT = ', np.mean(mean_EAT)

    print 'Differences over 50 mum:'
    print 'Delta PO2 max = ', mean_PO2_max[:-5] - mean_PO2_max[-4:]
    print 'Delta PO2 min = ', mean_PO2_min[:-5] - mean_PO2_min[-4:]
    print 'Delta PO2 mean = ', mean_PO2_mean[:-5] - mean_PO2_mean[-4:]

## Plot the PO2 pulse amplitude at various probes
#  @param minMaxDicts List with dictionaries read by extractProbeMinMax.
#  @param probes List of dictionaries produces by loadProbesEuler.py
def plotAveragedProbePulses(minMaxDicts, probes):
    
    nProbes = len(probes)
    
    # extract probe positions (scaled to mum)
    xProbes = [1e6*probe['positions'][0][0] for probe in probes]

    PO2_min  = [minMaxDict['PO2_min'] for minMaxDict in minMaxDicts]
    PO2_max  = [minMaxDict['PO2_max'] for minMaxDict in minMaxDicts]

    mean_PO2_max  = np.asarray([np.mean(P) for P in PO2_max])
    mean_PO2_min  = np.asarray([np.mean(P) for P in PO2_min])

    pulseAmplitude = mean_PO2_max - mean_PO2_min

    style = {'markerfacecolor': 'None',
             'markersize'     : 5,
             'linewidth'      : 0.3}

    plt.plot(xProbes, pulseAmplitude, 'gx-', label=r'pulse amplitude', **style)
    print 'Pulse amplitudes = ', pulseAmplitude


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
    figOptions = FigureOptions(parser)

    args            = parser.parse_args()
    fieldName       = args.field
    probeNames      = args.probeNames
    allProbes       = args.allProbes
    from_time       = args.from_time
    nearWallPulses  = args.nearWallPulses
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    if allProbes:
        probeNames = probeUtils.probeNames('.', suffix='PO2')

    probes = [loadProbes('.', probeName, fieldName) for probeName in probeNames]
    EAT_dicts = [extractEATs('.', probe, 0) for probe in probes]
    
    plotAveragedProbes(EAT_dicts, probes, from_time)

    if nearWallPulses:
        minMaxDicts = [extractProbeMinMax('.', probe, probeIdx=2,
                                          from_time=from_time) for probe in probes]
        plotAveragedProbePulses(minMaxDicts, probes)

    figOptions.adjustAxes()
    figOptions.setGrid()

    plotName = 'compare%iProbes' % len(probeNames)
    figOptions.saveFig(plotName)


