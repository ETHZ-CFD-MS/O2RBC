#!/usr/bin/env python
"""
Plot data extracted from several parameter studies for diffusive interaction between capillaries.
This script is used in the COSH paper.
"""

import argparse
import os

import matplotlib.pyplot as plt

from HbO2.plot.diffusiveinteraction import DiffusiveInteractionParameterStudyPlotter
from HbO2.plot.styles import set_COSH_rc_params
from HbO2.postprocess.factory.case import load_settings_file
from HbO2.postprocess.parameterstudy import ParameterStudyPostProcessor
from HbO2.setup.parameterStudy import ParameterStudyFactory
from plot.figureoptions import FigureOptions
from plot.utils import annotate_axis_corner

set_COSH_rc_params()

path_to_parameter_studies = "/local/aluecker/OpenFOAM/aluecker-2.3.0/run/cbf/HbO2/eulerGraph/COSH/straightCapillaryArray/diffusiveInteraction"

studies = ["domainHeightStudy",
           "M0Study",
           "RBCVelocityStudy",
           "LDMeanStudy"]
annotations = ['A', 'B', 'C', 'D']

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

fig, axs = plt.subplots(2, 2, sharey=True)
flattened_axs = [ax for sublist in axs for ax in sublist]

for study, ax, text in zip(studies, flattened_axs, annotations):
    path_to_study = os.path.join(path_to_parameter_studies, study)
    path_to_study_param_file = os.path.join(path_to_study, file_name)
    param_study = ParameterStudyFactory().makeParameterStudy(path_to_study_param_file)
    settings_dict = load_settings_file(path_to_study)
    postprocessor = ParameterStudyPostProcessor(param_study, settings_dict)
    plotter = DiffusiveInteractionParameterStudyPlotter(postprocessor, fig_options, model_names)
    plt.sca(ax)
    getattr(plotter, method_name)()
    annotate_axis_corner(ax, text)

axs[0,1].set_ylabel('')
axs[1,1].set_ylabel('')

fig_options.saveFig('plotDiffusiveInteractionMulti_{:s}'.format(method_name))
