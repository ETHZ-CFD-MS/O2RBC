#!/usr/bin/env python
"""
Plot longitudinal PO2 profiles in the tissue using the ODE model for several
cases.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric
from plot.figureoptions import FigureOptions
from plot.plotPO2ProfileAnalytical import plotPO2Profile


parser = argparse.ArgumentParser()
parser.add_argument('--cases', nargs='+', help='Relative paths of folders to plot',
                    default=['.'])
parser.add_argument('-r', type=float, help='Radius for PO2 calculation',
                    default = 17.6e-6)
fig_options = FigureOptions(parser)
args = parser.parse_args()
fig_options.parseOptions(args)
cases = args.cases
radius = args.r

styles = [{'color': 'k'},
          {'color': 'b'},
          {'color': 'r'}]

for case, style in zip(cases, styles):
    sim_params = IOHbO2ParametersAxisymmetric(case)
    plotPO2Profile(sim_params, radius=radius, style=style)

plt.legend([case.replace('_', '-') for case in cases])

plot_name = 'plot_compare_' + '_'.join(cases) + '_r_{:g}'.format(1e6*radius)
fig_options.saveFig(plot_name)

