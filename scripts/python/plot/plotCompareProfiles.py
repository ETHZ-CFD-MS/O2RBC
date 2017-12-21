#!/usr/bin/env python
#
# Plot on the same graph: profiles from different simulations
#

import os
import numpy as np
import argparse
import matplotlib.pyplot as plt

from plot.figureoptions import FigureOptions
from plotProfile import plotProfile

def setLegendCapArtDilationStudy():
    legends = ['baseline', 'cap. dilation', 'CBF+50\%', 'CBF \& cap. dil.']
    ax = plt.gca()
    box = ax.get_position()
    # for legend outside the plot, on the right
    # ax.set_position([box.x0, box.y0, box.width * 0.78, box.height])
    # ax.legend(legends, loc='upper center', fontsize=8, bbox_to_anchor=(1.17, 0.7), 
              # ncol=1, handletextpad=0.5, handlelength=2.4, labelspacing=1, fancybox=True, shadow=True)
    # for legend outside the plot, on the top
    ax.set_position([box.x0, box.y0, box.width * 0.78, box.height])
    ax.legend(legends, loc='upper center', fontsize=8, bbox_to_anchor=(1.17, 0.7), 
              ncol=1, handletextpad=0.5, handlelength=2.4, labelspacing=1, fancybox=True, shadow=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', '-p', help='Name of the profile file')
    parser.add_argument('--time', '-t', help='Time for which to plot the profile')
    parser.add_argument('--cases', help='Name of the cases to compare', 
                        nargs='+')
    figOptions = FigureOptions(parser)

    args  = parser.parse_args()
    profileName = args.profile
    time        = args.time
    cases       = args.cases
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    styles = []
    styles.append({'color': 'b',
                  'linestyle': '-'})
    styles.append({'color': 'r',
                  'linestyle': '--',
                  'dashes': (4,2.5)})
    styles.append({'color': 'k',
                  'linestyle': ':',
                  'dashes': (1,2)})
    styles.append({'color': 'g',
                  'linestyle': '-.',
                  'dashes': (4,2,1,2)})

    for i, case in enumerate(cases):
        domainPath = os.path.join(case, 'domain')
        plotProfile(domainPath, time, profileName, styles[i])
    # plt.ylim([18, 47])
    plt.ylim([15, 50])
    figOptions.adjustAxes()
    # setLegendCapArtDilationStudy()

    plotName = 'compareProfiles_%s' % profileName
    stringList = cases
    stringList.insert(0, plotName)
    plotName = '-'.join(stringList)
    figOptions.saveFig(plotName)
