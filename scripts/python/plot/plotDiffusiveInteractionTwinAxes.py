#!/usr/bin/env python
"""
Plot data extracted from several parameter studies for diffusive interaction between capillaries.
This script is used in the COSH paper.
"""

import argparse
import matplotlib.pyplot as plt
import os
import string

from HbO2.setup.parameterStudy import ParameterStudyFactory
from HbO2.plot import styles
from HbO2.plot.diffusiveinteraction import DiffusiveInteractionParameterStudyPlotter
from HbO2.postprocess.diffusiveinteraction import DiffusiveInteractionParamStudyPostProcessor
from plot.figureoptions import FigureOptions
from plot.utils import annotate_axis_corner

path_to_parameter_studies = "/local/aluecker/OpenFOAM/aluecker-2.3.0/run/cbf/HbO2/eulerGraph/COSH/straightCapillaryArray/diffusiveInteraction"

studies = ["RBCVelocityStudy",
           "LDMeanStudy",
           "M0Study",
           "domainHeightStudy"]
annotations = string.ascii_uppercase

parser = argparse.ArgumentParser()
parser.add_argument('--paramFile',
                    help='Path to parameter study file', default='params.json')
parser.add_argument('--methodName', '-m', help='Name of the plotting method to call')
parser.add_argument('--modelNames', nargs='+', help='Models for which to plot results',
                    default=['krogh', 'simple'])
fig_options = FigureOptions(parser)
args = parser.parse_args()
file_name = args.paramFile
method_name = args.methodName
model_names = args.modelNames
fig_options.parseOptions(args)

fig, axs = plt.subplots(1, 2, sharey=True)
all_axs = [[ax, ax.twiny()] for ax in axs]
flattened_axs = [ax for sublist in all_axs for ax in sublist]
style_scheme = styles.StyleSchemeCOSH()
twin_axes_color = style_scheme['twinAxis']['color']
colors = ['k', twin_axes_color, 'k', twin_axes_color]

for study, ax, color in zip(studies, flattened_axs, colors):
    path_to_study_file = os.path.join(path_to_parameter_studies, study, file_name)
    param_study = ParameterStudyFactory().makeParameterStudy(path_to_study_file)
    postprocessor = DiffusiveInteractionParamStudyPostProcessor(param_study)
    plotter = DiffusiveInteractionParameterStudyPlotter(postprocessor, fig_options, model_names)
    plt.sca(ax)
    getattr(plotter, method_name)(color=color)
    ax.xaxis.label.set_color(color)
    for t1 in ax.get_xticklabels():
        t1.set_color(color)

axs[1].set_ylabel('')
annotate_axis_corner(axs[0], 'A')
annotate_axis_corner(axs[1], 'B')

fig_options.saveFig('plotDiffusiveInteractionTwinAxes_{:s}'.format(method_name))
