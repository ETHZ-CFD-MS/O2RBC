#! /usr/bin/env python
"""
Setup, postprocess and plot a named parameter study for diffusive interaction
"""

import argparse

from HbO2.COSH.plotCOSH import LDvRBCCOSHParamStudyPlotter
from HbO2.plot.styles import set_COSH_rc_params
from HbO2.postprocess.factory.parameterstudy import make_param_study_post_processor
from HbO2.setup.parameterStudy import ParameterStudyFactory
from plot.figureoptions import FigureOptions

set_COSH_rc_params()


parser = argparse.ArgumentParser()
parser.add_argument('--paramFile', help='Path to parameter study file', default='params.json')
parser.add_argument('--manual-labels', action='store_true',
                    help='Whether to set manually the contour labels')
fig_options = FigureOptions(parser)
args = parser.parse_args()
manual_labels = args.manual_labels
file_name = args.paramFile
fig_options.parseOptions(args)

param_study = ParameterStudyFactory().makeParameterStudy(file_name)
postprocessor = make_param_study_post_processor(param_study)
plotter = LDvRBCCOSHParamStudyPlotter(postprocessor, fig_options)
if fig_options.blackWhite:
    plotter.cmap = None
else:
    plotter.colors = None
plotter.plotMultiPanelKOSContourAndOscillationDistance(manual_labels)
