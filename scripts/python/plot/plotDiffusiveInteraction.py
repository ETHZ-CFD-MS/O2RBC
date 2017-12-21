#!/usr/bin/env python
"""
Plot Hb profiles in a capillary using the analytical model
based on Krogh's assumptions
"""

import argparse

import matplotlib.pyplot as plt

from HbO2.plot.diffusiveinteraction import DiffusiveInteractionComparisonPlotter, \
                                           DiffusiveInteractionModelPlotter
from HbO2.postprocess.factory.case import make_post_processor
from HbO2.setup.geometry import ParallelCapillaries
from HbO2.setup.simulationParameters import IOHbO2ParametersStraightCapillaries
from HbO2.plot.styles import create_COSH_legend
from plot.figureoptions import FigureOptions
from plot.utils import add_multi_panel_parser_options, annotate_axis_corner


parser = argparse.ArgumentParser()
parser.add_argument('--modelOnly', '-m', help='Plot only model results',
                    action='store_true')
parser.add_argument('--modelNames', nargs='+', help='Models for which to plot results',
                    default=['krogh', 'simple', 'equal_outfluxes', 'linearized_ODE'])
parser.add_argument('--exponentialFit', help='Whether to plot the exponential fit',
                    action='store_true')
parser.add_argument('--noMean', help='Whether to plot the mean profiles',
                    action='store_false')
add_multi_panel_parser_options(parser)
figOptions = FigureOptions(parser)
args = parser.parse_args()
model_only = args.modelOnly
model_names = args.modelNames
exponential_fit = args.exponentialFit
plot_mean = not args.noMean
linearized_ode = 'linearized_ODE' in model_names
try:
    model_names.remove('linearized_ODE')
except ValueError:
    pass
figOptions.parseOptions(args)

parallel_capillaries = ParallelCapillaries.fromKeyValueFile('geometricData', 2)
simParams = IOHbO2ParametersStraightCapillaries('.')
postprocessor = make_post_processor('.')

plot_name_suffix = '' if simParams.cocurrentFlow() else '_countercurrent'
if model_only:
    for model_name in postprocessor.integrators:
        plotter = DiffusiveInteractionModelPlotter(postprocessor.integrators[model_name])
        plotter.plotAll()
        figOptions.saveFig('plot_{:s}_modelDiffusiveInteraction_Hb{:s}'.\
                           format(model_name, plot_name_suffix))
        plt.clf()
else:
    plotter = DiffusiveInteractionComparisonPlotter(postprocessor, model_names)
    if args.multipanel:
        panel_layout = (2, 1) if not args.panelLayout else args.panelLayout
        fig, axs = plt.subplots(*panel_layout, sharex=True)
        plt.sca(axs[0])
        plotter.plotHbProfiles(plot_mean=plot_mean)
        annotate_axis_corner(axs[0], 'A')
        if panel_layout[1] == 1:
            axs[0].set_xlabel('')
        plt.sca(axs[1])
        plotter.plotHbDiffProfiles(exponential_fit=exponential_fit, linearized_ode=linearized_ode)
        # legend = create_COSH_legend(axs[1], bbox_to_anchor=(1, 0.86))
        legend = create_COSH_legend(axs[1], bbox_to_anchor=(0.5, 0.86))
        annotate_axis_corner(axs[1], 'B')
        figOptions.saveFig('plotDiffusiveInteractionMultipanel' + plot_name_suffix)
    else:
        plotter.plotHbProfiles()
        plotter.plotHbDiffProfiles(exponential_fit=exponential_fit)
        figOptions.saveFig('plotDiffusiveInteraction' + plot_name_suffix)
        plt.clf()
        for model_name, modelPlotter in plotter.modelPlotters.iteritems():
            modelPlotter.plotOutfluxes()
        figOptions.saveFig('plotOutfluxes')
        plt.clf()

        for model_name, modelPlotter in plotter.modelPlotters.iteritems():
            modelPlotter.plotTissueRadii()
        figOptions.saveFig('plotTissueRadii')
        plt.clf()

        plotter.plotDistalHbTemporalProfile()
        figOptions.saveFig('plotDistalHb')
