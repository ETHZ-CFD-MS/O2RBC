#!/usr/bin/env python
#
# Plot simulated and analytical longitudinal profiles of Hb.
#

import argparse

from HbO2.plot.simulated import AxisymmetricCasePlotter, GraphCasePlotter
from HbO2.postprocess.case import GraphCasePostProcessor
from HbO2.postprocess.factory.case import make_post_processor
from HbO2.setup.case import SimulationParametersFactory
from HbO2.setup.utils import isAxisymmetricCase, isGraphCase
from plot.figureoptions import FigureOptions
from plotHbProfileAnalytical import plotHbProfile


def plotHbSimulatedAxisymmetric(case_name):
    postprocessor = make_post_processor(case_name)
    plotter = AxisymmetricCasePlotter(postprocessor)
    plotter.plot_field_average_profile('Hb_mean', n_average=1)


def plotHbSimulatedGraph(caseName):
    postprocessor = GraphCasePostProcessor(caseName)
    plotter = GraphCasePlotter(postprocessor)
    plotter.plotFieldAverageProfile('Hb_mean', 0, nAverage=10)


def plotHbSimulatedAndAnalytical(caseName, simParams, **kwargs):
    if isGraphCase(caseName):
        plotHbSimulatedGraph(caseName)
    elif isAxisymmetricCase(caseName):
        plotHbSimulatedAxisymmetric(caseName)
    plotHbProfile(simParams, **kwargs)
    # plt.legend(['Simulated', 'ODE model'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    figOptions = FigureOptions(parser)
    args = parser.parse_args()
    figOptions.parseOptions(args)

    analyticalLineStyle = {'color': 'k',
                           'linestyle': '-',
                           'dashes': (4, 2.5)}

    simParams = SimulationParametersFactory().make_sim_params('.')
    plotHbSimulatedAndAnalytical('.', simParams, style=analyticalLineStyle)
    figOptions.saveFig('plotHbSimulatedAndAnalytical')
