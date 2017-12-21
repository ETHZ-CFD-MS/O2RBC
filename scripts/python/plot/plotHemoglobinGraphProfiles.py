#!/usr/bin/env python
"""
Plot hemoglobin saturation profiles for one OpenFOAM simulation in a capillary network.
"""

import argparse
import matplotlib.pyplot as plt
import string

from utilities.arguments import add_case_argument

from HbO2.postprocess.factory.case import make_post_processor
from plot.utils import annotate_axis_corner

from HbO2.plot.hemoglobingraph import HemoglobinOnSegmentsPlotter
from plot.figureoptions import FigureOptions

annotations = string.ascii_uppercase
ncol = 2

parser = argparse.ArgumentParser()
add_case_argument(parser)
parser.add_argument('--settingsFile', help='Relative path to the file with postprocessing settings',
                    default='postprocess.json')
edge_group = parser.add_mutually_exclusive_group()
edge_group.add_argument('--allSegments', '-a', action='store_true',
                        help='Plot all segments that have at least nProfiles paths')
edge_group.add_argument('--segmentIndices', '-s', type=int, nargs='+', default=[],
                        help='Segment indices where to plot hemoglobin saturation')
parser.add_argument('--integratedProfiles', action='store_true')
parser.add_argument('--nProfileMin', '-n', type=int, default=2,
                    help='Minimal number of profiles to plot')
parser.add_argument('--multipanel', action='store_true',
                    help='Whether to create multiple panels in one figure')
figOptions = FigureOptions(parser)
args = parser.parse_args()
case_name = args.caseName
settings_file = args.settingsFile
integrated_profiles = args.integratedProfiles
n_profile_min = args.nProfileMin
multipanel = args.multipanel
figOptions.parseOptions(args)

postprocessor = make_post_processor(case_name, param_file=settings_file)
plotter = HemoglobinOnSegmentsPlotter(postprocessor, n_profile_min=n_profile_min)

if args.allSegments:
    sids = plotter.edge_ids()
else:
    sids = args.segmentIndices

if multipanel:
    nrow = len(sids)/ncol
    fig, axs = plt.subplots(nrow, ncol)
    flattened_axs = [ax for sublist in axs for ax in sublist]
    for si, ax, annotation in zip(sids, flattened_axs, annotations):
        plt.sca(ax)
        plotter.plot_hb_profiles(si)
        plotter.plot_hb_mean_pm_std(si)
        plotter.restrict_x_limits_to_defined_values(si)
        annotate_axis_corner(ax, annotation)
    e_string = '_'.join(['{:d}'.format(si) for si in sids])
    plotname = 'plotHbSimulatedEdges_e_{:s}'.format(e_string)
    figOptions.saveFig(plotname)
else:
    for si in sids:
        plt.clf()
        plotter.plot_hb_profiles(si)
        plotter.plot_hb_mean_pm_std(si)
        if integrated_profiles:
            plotter.plot_hb_profile_integrated(si)
        plotter.restrict_x_limits_to_defined_values(si)
        n_profiles = plotter.n_profiles(si)
        plotname = 'plotHbSimulatedEdges_e_{:d}_n_{:d}'.format(si, n_profiles)
        figOptions.saveFig(plotname)
        plt.clf()
        plotter.plot_hb_drop_vs_time_to_previous_rbc(si)
        plotname = 'plotHbDropVsTimeToPreviousRBC_e_{:d}'.format(si, n_profiles)
        figOptions.saveFig(plotname)
        plt.clf()
        s = postprocessor.rbcDataPostProcessor.rbc_path_analyzer.s_coord_with_most_paths(si)
        n_profiles = postprocessor.n_rbc_average(si)
        plotter.plotFieldTemporalProfile('Hb_mean', s, si, n_profiles)
        plotname = 'plotHbTemporal_e_{:d}'.format(si)
        figOptions.saveFig(plotname)