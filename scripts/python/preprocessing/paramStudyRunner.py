#! /usr/bin/env python
"""
Setup, postprocess and plot a named parameter study for diffusive interaction
"""

import argparse


from HbO2.COSH.postprocess import COSHParamStudyPostProcessor,\
                                  COSHLDvRBCParamStudyPostProcessor
from HbO2.plot.factory.parameterstudy import make_param_study_plotter
from HbO2.postprocess.factory.case import load_settings_file
from HbO2.postprocess.factory.parameterstudy import make_param_study_post_processor
from HbO2.postprocess.parameterstudy import LDvRBCParameterStudyPostProcessor,\
                                            ParameterStudyPostProcessorWriter
from HbO2.setup.parameterStudy import ParameterStudyFactory
from plot.figureoptions import FigureOptions


parser = argparse.ArgumentParser()
parser.add_argument('--setup', '-s', action='store_true')
parser.add_argument('--postProcess', '-p', action='store_true')
parser.add_argument('--plot', action='store_true')
parser.add_argument('--paramFile',
    help='Path to parameter study file', default='params.json')
parser.add_argument('--settingsFile', help='Relative path to the file with postprocessing settings',
                    default='postprocess.json')
parser.add_argument('--counter', '-c', help='Countercurrent flow',
    action='store_true')
fig_options = FigureOptions(parser)
args = parser.parse_args()
doSetup = args.setup
doPostProcess = args.postProcess
doPlot = args.plot
fileName = args.paramFile
settings_file = args.settingsFile
cocurrent = not args.counter

param_study = ParameterStudyFactory().makeParameterStudy(fileName)

if doSetup:
    param_study.setup()
else:
    post_processor = make_param_study_post_processor(param_study)
    if doPostProcess:
        post_processor.run()
        writer = ParameterStudyPostProcessorWriter(post_processor)
        writer.write_results()
    elif doPlot:
        fig_options.parseOptions(args)
        plotter = make_param_study_plotter(post_processor, fig_options)
        plotter.plot_all()
