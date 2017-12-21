#!/usr/bin/env python
"""
Plot statistics of functional and topological volumes in a capillary network.
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np

from HbO2.plot import styles as cosh_styles
from HbO2.plot.tissuevolumes import TissueVolumePlotter
from HbO2.postprocess.results import MultiCaseResults
from HbO2.postprocess.tissuevolumes import TissueVolumeResults

from plot.figureoptions import FigureOptions
from utilities.arguments import add_single_or_multi_case_group

style_scheme = cosh_styles.StyleSchemeCOSH()
plot_styles = [style_scheme['MVN1'], style_scheme['MVN2']]

parser = argparse.ArgumentParser()
add_single_or_multi_case_group(parser)
fig_options = FigureOptions(parser)
args = parser.parse_args()
multi_case = True if args.caseNames else False
if not multi_case:
    case_paths = [args.caseName]
else:
    case_paths = args.caseNames
fig_options.parseOptions(args)

results = [TissueVolumeResults(case_path) for case_path in case_paths]
plotters = [TissueVolumePlotter(res) for res in results]
combined_results = MultiCaseResults(results, np.hstack)
combined_plotter = TissueVolumePlotter(combined_results)
styles = [dict() for _ in plotters]
labels = []
for case_path, style in zip(case_paths, styles):
    if '20151213' in case_path or 'MVN1' in case_path:
        style.update(style_scheme['MVN1'])
        style.update({'label': 'MVN 1'})
    elif '20150715' in case_path or 'MVN2' in case_path:
        style.update(style_scheme['MVN2'])
        style.update({'label': 'MVN 2'})

plt.clf()
combined_plotter.plot_histogram_tissue_radii()
fig_options.saveFig('plotHistogramTissueRadii')

plt.clf()
for plotter, style in zip(plotters, styles):
    plotter.plot_functional_vs_topological_tissue_radii(**style)
cosh_styles.create_COSH_legend(plt.gca(), loc='upper left', markerscale=0.7, handletextpad=0,
                               handlelength=2.5)
fig_options.saveFig('plotFunctionalVsTopologicalTissueRadii')

plt.clf()
for plotter, style in zip(plotters, styles):
    plotter.plot_tissue_radii_vs_hb_mean(**style)
cosh_styles.create_COSH_legend(plt.gca(), loc='upper left', markerscale=0.8, handletextpad=0,
                               handlelength=2.5)
fig_options.saveFig('plotTissueRadiiVsHbMean')

plt.clf()
for plotter, style in zip(plotters, styles):
    plotter.plot_tissue_radii_vs_hb_drop(**style)
fig_options.saveFig('plotTissueRadiiVsHbDrop')

plt.clf()
for plotter, style in zip(plotters, styles):
    plotter.plot_tissue_radii_vs_rbc_flow(**style)
cosh_styles.create_COSH_legend(plt.gca(), loc='best', markerscale=0.8, handletextpad=0,
                               handlelength=2.5)
fig_options.saveFig('plotTissueRadiiVsRBCFlow')

plt.clf()
for plotter, style in zip(plotters, styles):
    plotter.plot_tissue_radii_vs_rbc_velocity(**style)
fig_options.saveFig('plotTissueRadiiVsRBCVelocity')

plt.clf()
for plotter, style in zip(plotters, styles):
    plotter.plot_tissue_radii_vs_hematocrit(**style)
fig_options.saveFig('plotTissueRadiiVsHematocrit')

plt.clf()
for plotter, style in zip(plotters, styles):
    plotter.plot_tissue_radii_vs_arrival_transit_time(**style)
fig_options.saveFig('plotTissueRadiiVsArrivalTransitTime')

plt.clf()
for plotter, style in zip(plotters, styles):
    plotter.plot_tissue_radii_vs_arrival_path_length(**style)
fig_options.saveFig('plotTissueRadiiVsArrivalPathLength')

plt.clf()
for plotter, style in zip(plotters, styles):
    plotter.plot_tissue_radius_difference_vs_hb_mean(**style)
fig_options.saveFig('plotTissueRadiusDifferenceVsHbMean')
