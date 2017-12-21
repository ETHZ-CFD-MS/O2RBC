#!/usr/bin/env python
"""
Plot Hb profiles in a capillary using the analytical model
based on Krogh's assumptions
"""

import argparse

import matplotlib.pyplot as plt

from HbO2.plot.simulated import GraphCasePlotter
from plot.figureoptions import FigureOptions

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    edge_group = parser.add_mutually_exclusive_group()
    edge_group.add_argument('--allEdges', '-a', action='store_true',
                            help='Plot all edges that have at least nProfiles paths')
    edge_group.add_argument('--edgeIndices', '-e', type=int, nargs='+',
                            help='Edge indices where to plot hemoglobin saturation')
    parser.add_argument('--nProfiles', '-n', type=int, default=100,
                        help='Maximal number of profiles to plot')
    figOptions = FigureOptions(parser)
    args = parser.parse_args()
    n_profiles = args.nProfiles
    figOptions.parseOptions(args)

    plotter = GraphCasePlotter('.')

    if args.allEdges:
        eids = plotter.postprocessor.rbcDataPostProcessor.edgeIdsWithPaths(nPathMin=n_profiles)
    else:
        eids = args.edgeIndices

    for ei in eids:
        plt.clf()
        plotter.plotFieldProfiles('Hb_mean', ei, n_profiles,
                                  style={'color': 'k', 'linewidth': 0.2})
        plotter.plotFieldAverageProfile('Hb_mean', ei, nAverage=n_profiles,
                                        style={'color': 'b'})
        plotter.plotFieldAverageWithStdProfile('Hb_mean', ei, nAverage=n_profiles,
                                               meanStyle={'color': 'k', 'linewidth': 2},
                                               stdStyle={'color': 'b', 'linewidth': 1})
        plotname = 'plotHbSimulatedEdges_e_{:d}_n_{:d}'.format(ei, n_profiles)
        figOptions.saveFig(plotname)
