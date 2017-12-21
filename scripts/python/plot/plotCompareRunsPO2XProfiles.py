#!/usr/bin/env python
#
# Plot on the same graph: x-profiles (snapshot) from different
# simulations
#

import numpy as np
import argparse

from plot.figureoptions import FigureOptions
from plotPO2XProfiles import plotPO2ProfileEulerian

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--time', type=float, help='Time for the snapshot')
    parser.add_argument('--cases', help='Name of the cases to compare', 
                        nargs='+')
    figOptions = FigureOptions(parser)

    args  = parser.parse_args()
    time  = args.time
    cases = args.cases
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    i = 0
    for case in cases:
        domainPath = '%s/domain' % case
        alpha = 1. - i/(float(len(cases)))
        plotPO2ProfileEulerian(domainPath, float(time), alpha)
        i = i+1

    figOptions.adjustAxes()
    figOptions.setGrid()
    plotName = 'comparePO2XProfiles_t_%.3f' % time
    stringList = cases
    stringList.insert(0, plotName)
    plotName = '-'.join(stringList)
    figOptions.saveFig(plotName)

