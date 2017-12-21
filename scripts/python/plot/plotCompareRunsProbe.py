#!/usr/bin/env python
#
# Plot on the same graph: probe profiles from different
# simulations
#

import argparse

import matplotlib.pyplot as plt

from plot.figureoptions import FigureOptions
from plotProbe import plot_probe

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--probeName', '-p', help='Name of the probe directory', 
                        default='probe05PO2')
    parser.add_argument('--positions', type=int, nargs='+', 
                        help='Positions to plot', default=[-1])
    parser.add_argument('--min-time', type=float, help='Smallest time that should be used',
                        default = 0.0)
    parser.add_argument('--max-time', type=float, help='Largest time that should be used',
                        default = 1.0)
    parser.add_argument('--cases', help='Name of the cases to compare', 
                        nargs='+')
    figOptions = FigureOptions(parser)

    args  = parser.parse_args()
    probeName = args.probeName
    probeIdx  = args.positions
    min_time  = args.min_time
    max_time  = args.max_time
    cases     = args.cases
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    fieldName = 'PO2'

    i = 0
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
    for case in cases:
        domainPath = '%s/domain' % case

        alpha = 1. - i/(float(len(cases)))
        plot_probe(domainPath, fieldName, probeName, probeIdx,
                   min_time, max_time, styles[i])
        i = i+1

    # plt.legend(['std LD = 0.0', 
                # 'std LD = 0.08',  
                # 'std LD = 0.16'])
    plt.ylim([0, 13])
    figOptions.adjustAxes()
    figOptions.setGrid()
    plotName = 'compare%s' % probeName
    stringList = cases
    stringList.insert(0, plotName)
    plotName = '-'.join(stringList)
    figOptions.saveFig(plotName)


