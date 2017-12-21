#!/usr/bin/env python
"""
Plot the hemoglobin saturation distribution on distal side for one OpenFOAM simulation
in a capillary network.
"""

import argparse

from HbO2.plot.styles import StyleSchemeCOSH
from HbO2.plot.hemoglobingraph import HemoglobinOnWholePathsPlotter
from HbO2.postprocess.factory.case import make_post_processor
from plot.figureoptions import FigureOptions
from utilities.arguments import add_case_argument

parser = argparse.ArgumentParser()
add_case_argument(parser)
parser.add_argument('--settingsFile', help='Relative path to the file with postprocessing settings',
                    default='postprocess.json')
parser.add_argument('--noErrorBars', action='store_true')
figOptions = FigureOptions(parser)
args = parser.parse_args()
case_name = args.caseName
settings_file = args.settingsFile
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

style_scheme = StyleSchemeCOSH()
sim_linestyle = style_scheme['individualRBCWholePath']
symbol_style = style_scheme['individualRBCWholePath']
int_linestyle = {'color': '0.2', 'alpha': 1.0}

annotations = ('A', 'B', 'C', 'D')
# annotations = ('E', 'F', 'G', 'H')
plotter.plot_compared_hb_distal_distribution_multipanel(panel_annotation=annotations)
figOptions.saveFig('plotComparedHbDistalDistributionMultiPanel' + suffix + average_suffix)
