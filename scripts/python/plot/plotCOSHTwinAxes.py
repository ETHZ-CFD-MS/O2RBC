#!/usr/bin/env python
"""
Plot data extracted from several parameter studies for COSH.
This script is used in the COSH paper.
"""

import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os
import string

from HbO2.COSH.plotCOSH import COSHParamStudyPlotter
from HbO2.plot import styles
from HbO2.postprocess.factory.parameterstudy import make_param_study_post_processor
from HbO2.setup.parameterStudy import ParameterStudyFactory
from plot.figureoptions import FigureOptions
from plot.utils import annotate_axis_corner, add_panel_layout_parser_options

styles.set_COSH_rc_params()


path_to_parameter_studies = "/local/aluecker/OpenFOAM/aluecker-2.3.0/run/cbf/HbO2/eulerAxisymmetric/COSH"

studies = ["LD_U_study/paramStudyShort2_slices/paramStudyShort2_LD_0.3",
           "LD_U_study/paramStudyShort2_slices/paramStudyShort2_U_0.001",
           "sigmaSStudy",
           "O2ConsumptionStudyShort"]
annotations = string.ascii_uppercase

parser = argparse.ArgumentParser()
parser.add_argument('--paramFile',
                    help='Path to parameter study file', default='params.json')
parser.add_argument('--methodName', '-m', help='Name of the plotting method to call')
add_panel_layout_parser_options(parser)
fig_options = FigureOptions(parser)
args = parser.parse_args()
file_name = args.paramFile
method_name = args.methodName
fig_options.parseOptions(args)

panel_layout = (2, 1) if not args.panelLayout else args.panelLayout
fig, axs = plt.subplots(*panel_layout, sharey=True)
all_axs = [[ax, ax.twiny()] for ax in axs]
flattened_axs = [ax for sublist in all_axs for ax in sublist]
style_scheme = styles.StyleSchemeCOSH()
twin_axes_color = style_scheme['twinAxis']['color']
colors = ['k', twin_axes_color, 'k', twin_axes_color]

for study, ax, color in zip(studies, flattened_axs, colors):
    path_to_study_file = os.path.join(path_to_parameter_studies, study, file_name)
    param_study = ParameterStudyFactory().makeParameterStudy(path_to_study_file)
    postprocessor = make_param_study_post_processor(param_study)
    plotter = COSHParamStudyPlotter(postprocessor, fig_options)
    plt.sca(ax)
    getattr(plotter, method_name)(style={'color': color})
    ax.xaxis.label.set_color(color)
    for t1 in ax.get_xticklabels():
        t1.set_color(color)

axs[0].xaxis.set_major_locator(MultipleLocator(0.4))
axs[0].yaxis.set_major_locator(MultipleLocator(10))
flattened_axs[1].xaxis.set_major_locator(MultipleLocator(0.1))
flattened_axs[2].xaxis.set_major_locator(MultipleLocator(0.02))
flattened_axs[3].xaxis.set_major_locator(MultipleLocator(0.2))

if panel_layout[0] == 1:
    axs[1].set_ylabel('')
annotate_axis_corner(axs[0], 'A')
annotate_axis_corner(axs[1], 'B')

fig_options.saveFig('plotCOSHTwinAxes_{:s}'.format(method_name))
