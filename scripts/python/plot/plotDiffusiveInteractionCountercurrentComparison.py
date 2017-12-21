#!/usr/bin/env python
"""
Plot data extracted two parameter studies for the same parameter, one with cocurrent flow
and the other with countercurrent flow.
"""

import argparse
import matplotlib.pyplot as plt
import os

from HbO2.postprocess.factory.case import load_settings_file, make_post_processor
from HbO2.postprocess.parameterstudy import ParameterStudyPostProcessor
from HbO2.plot.diffusiveinteraction import DiffusiveInteractionComparisonPlotter,\
                                           DiffusiveInteractionParameterStudyPlotter
from HbO2.setup.parameterStudy import ParameterStudyFactory
from plot.figureoptions import FigureOptions
from plot.utils import add_multi_panel_parser_options, annotate_axis_corner


def add_arrows(ax):
    ax.annotate("",
                xytext = (40, 0.735),
                textcoords='data',
                xy = (60, 0.705),
                xycoords='data',
                arrowprops=dict(arrowstyle='->', color='k',
                connectionstyle="arc3"),
    )
    ax.annotate("",
                xytext = (60, 0.58),
                textcoords='data',
                xy = (40, 0.57),
                xycoords='data',
                arrowprops=dict(arrowstyle='->', color='k',
                                connectionstyle="arc3"),
    )


parameter_studies = [("/local/aluecker/OpenFOAM/aluecker-2.3.0/run/cbf/HbO2/eulerGraph/COSH"
                      "/straightCapillaryArray/diffusiveInteraction/deltaSStudy"),
                     ("/local/aluecker/OpenFOAM/aluecker-2.3.0/run/cbf/HbO2/eulerGraph/COSH"
                      "/straightCapillaryArray/diffusiveInteraction/countercurrentFlow/deltaSStudy")]
case_name_profile = "case_deltaS_0.2"
model_names = ['krogh', 'equal_outfluxes']

parser = argparse.ArgumentParser()
parser.add_argument('--paramFile',
                    help='Path to parameter study file', default='params.json')
parser.add_argument('--methodName', '-m', help='Name of the plotting method to call')
add_multi_panel_parser_options(parser)
fig_options = FigureOptions(parser)
args = parser.parse_args()
file_name = args.paramFile
method_name = args.methodName
fig_options.parseOptions(args)

case_names = [os.path.join(param_study_path, case_name_profile)
              for param_study_path in parameter_studies]

panel_layout = (2, 1) if not args.panelLayout else args.panelLayout
fig, axs = plt.subplots(*panel_layout)

plt.sca(axs[0])
case_name = os.path.join(parameter_studies[1], case_name_profile)
postprocessor = make_post_processor(case_name)
plotter = DiffusiveInteractionComparisonPlotter(postprocessor, model_names)
plotter.plotHbProfiles()
add_arrows(axs[0])

plt.sca(axs[1])
styles = [{'style': {'linestyle': '-', 'color': 'k'}},
          {'style': {'linestyle': '--', 'color': 'k', 'dashes': (4, 4)}}]
for param_study_path, style in zip(parameter_studies, styles):
    param_study_file = os.path.join(param_study_path, file_name)
    param_study = ParameterStudyFactory().makeParameterStudy(param_study_file)
    settings_dict = load_settings_file(param_study_path)
    postprocessor = ParameterStudyPostProcessor(param_study, settings_dict)
    plotter = DiffusiveInteractionParameterStudyPlotter(postprocessor, fig_options, model_names)
    getattr(plotter, method_name)(**style)

annotate_axis_corner(axs[0], 'A')
annotate_axis_corner(axs[1], 'B')

fig_options.saveFig('plotDiffusiveCountercurrentComparison_{:s}'.format(method_name))
