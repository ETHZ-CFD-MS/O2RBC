#!/usr/bin/env python
"""
Plot hemoglobin saturation profiles for one OpenFOAM simulation in a capillary network.
"""

import argparse
import matplotlib.pyplot as plt
import string

from HbO2.plot.styles import StyleSchemeCOSH
from HbO2.plot.hemoglobingraph import HemoglobinOnWholePathsPlotter
from HbO2.postprocess.factory.case import make_post_processor
from plot.figureoptions import FigureOptions
from utilities.arguments import add_case_argument

annotations = string.ascii_uppercase

parser = argparse.ArgumentParser()
add_case_argument(parser)
parser.add_argument('--settingsFile', help='Relative path to the file with postprocessing settings',
                    default='postprocess.json')
parser.add_argument('--integratedProfiles', action='store_true')
parser.add_argument('--noErrorBars', action='store_true')
figOptions = FigureOptions(parser)
args = parser.parse_args()
case_name = args.caseName
settings_file = args.settingsFile
integrated_profiles = args.integratedProfiles
error_bars = not args.noErrorBars
figOptions.parseOptions(args)

postprocessor = make_post_processor(case_name, param_file=settings_file)
plotter = HemoglobinOnWholePathsPlotter(postprocessor)

path_ids = postprocessor.selected_path_indices()

if postprocessor.selection_mode == 'nPath':
    suffix = '_n_{:d}'.format(postprocessor.n_path)
elif postprocessor.selection_mode == 'firstTime':
    suffix = '_fromTime_{:g}'.format(postprocessor.selection_first_time)
else:
    suffix = ''
average_suffix = '_hbAveraged' if postprocessor.integrator.average_hb else ''
integrated_suffix = '_withIntegrated' if integrated_profiles else ''

style_scheme = StyleSchemeCOSH()
sim_linestyle = style_scheme['individualRBCWholePath']
symbol_style = style_scheme['individualRBCWholePath']
int_linestyle = {'color': '0.2', 'alpha': 1.0}

if integrated_profiles:
    plotter.plot_integrated_hb_profile_vs_path_length(style=int_linestyle)
plotter.plot_distal_hb_vs_path_length_with_paths(path_ids, path_style=sim_linestyle,
                                                 symbol_style=style_scheme['distalHbMarker'],
                                                 error_bars=error_bars)
figOptions.saveFig('plotHbWholePathVsPathLength' + suffix + integrated_suffix)

plt.clf()
if integrated_profiles:
    plotter.plot_integrated_hb_profile_vs_transit_time(style=int_linestyle)
plotter.plot_distal_hb_vs_transit_time_with_paths(path_ids, path_style=sim_linestyle,
                                                  symbol_style=style_scheme['distalHbMarker'],
                                                  error_bars=error_bars)
figOptions.saveFig('plotHbWholePathVsTransitTime' + suffix + integrated_suffix)

plt.clf()
plotter.plot_hb_drop_vs_transit_time()
figOptions.saveFig('plotHbDropVsTransitTime' + suffix)

plt.clf()
plotter.plot_hb_drop_vs_proximal_hb()
figOptions.saveFig('plotHbDropVsProximalHb' + suffix)

plt.clf()
plotter.plot_mean_hb_slope_transit_time_vs_proximal_hb()
figOptions.saveFig('plotHbSlopeTransitTimeVsProximalHb' + suffix)

plt.clf()
plotter.plot_mean_hb_slope_path_length_vs_proximal_hb()
figOptions.saveFig('plotHbSlopePathLengthVsProximalHb' + suffix)

plt.clf()
plotter.plot_compared_hb_distal_distribution()
figOptions.saveFig('plotComparedHbDistalDistribution' + suffix + average_suffix)

plt.clf()
plotter.plot_compared_hb_distal_distribution_multipanel()
figOptions.saveFig('plotComparedHbDistalDistributionMultiPanel' + suffix + average_suffix)
