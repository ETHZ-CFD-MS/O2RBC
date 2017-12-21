#!/usr/bin/env python
"""
Plot data extracted from several parameter studies for COSH.
This script is used in the COSH paper.
"""

import argparse
import os

import matplotlib.pyplot as plt


from HbO2.COSH.plotCOSH import COSHParamStudyPlotter
from HbO2.plot.styles import set_COSH_rc_params
from HbO2.postprocess.factory.parameterstudy import make_param_study_post_processor
from HbO2.setup.parameterStudy import ParameterStudyFactory
from plot.figureoptions import FigureOptions
from plot.utils import annotate_axis_corner

set_COSH_rc_params()

path_to_parameter_studies = "/local/aluecker/OpenFOAM/aluecker-2.3.0/run/cbf/HbO2/eulerAxisymmetric/COSH"
studies = ["O2ConsumptionStudyShort",
           "sigmaSStudy",
           "LD_U_study/paramStudyShort2_slices/paramStudyShort2_LD_0.3",
           "LD_U_study/paramStudyShort2_slices/paramStudyShort2_U_0.001"]
annotations = ['A', 'B', 'C', 'D']

parser = argparse.ArgumentParser()
parser.add_argument('--paramFile',
                    help='Path to parameter study file', default='params.json')
parser.add_argument('--methodName', '-m', help='Name of the plotting method to call')
fig_options = FigureOptions(parser)
args = parser.parse_args()
file_name = args.paramFile
method_name = args.methodName
fig_options.parseOptions(args)

fig, axs = plt.subplots(2, 2, sharey=True)
flattened_axs = [ax for sublist in axs for ax in sublist]

for study, ax, text in zip(studies, flattened_axs, annotations):
    path_to_study_file = os.path.join(path_to_parameter_studies, study, file_name)
    param_study = ParameterStudyFactory().makeParameterStudy(path_to_study_file)
    postprocessor = make_param_study_post_processor(param_study)
    plotter = COSHParamStudyPlotter(postprocessor, fig_options)
    plt.sca(ax)
    getattr(plotter, method_name)()
    annotate_axis_corner(ax, text)

axs[0,1].set_ylabel('')
axs[1,1].set_ylabel('')

fig_options.saveFig('plotCOSHMulti_{:s}'.format(method_name))

