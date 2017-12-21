#!/usr/bin/env python
"""
Plot slopes of PO2 with comparison to experimental data
"""

import argparse
import matplotlib.pyplot as plt
from plot.figureoptions import FigureOptions
from HbO2.plot import labels

parser = argparse.ArgumentParser()
fig_options = FigureOptions(parser)
args = parser.parse_args()
fig_options.parseOptions(args)

markersize = 8

# The plotted differences were obtained with the command
# plotLongitudinalGradients.py -a -f PO2 --cases final_Rl_16.1_Rr_16.1 final_tol_1e-3
# in /local/aluecker/OpenFOAM/aluecker-2.1.1/run/cbf/HbO2/eulerAxisymmetric/Parpaleix2013_validation
plt.plot(0, 12.25, '.', color='k', markersize=markersize)
plt.plot(0, 5.2, '.', color='b', markersize=markersize)
plt.plot(0, 3.0, '.', color=(0.86, 0.5, 0), markersize=markersize)
plt.errorbar(0, 3.0, yerr=2.7, color=(0.86, 0.5, 0))
labels.setYLabel('\Delta \mathrm{PO}_2', '\mathrm{mm\,Hg}/(50\,\muup m)')
plt.gca().axes.get_xaxis().set_visible(False)

plot_name = 'PO2SlopeComparisonExperiments'
fig_options.saveFig(plot_name)
