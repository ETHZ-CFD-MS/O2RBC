#!/usr/bin/env python
"""
Plot hemoglobin saturation profiles and related quantities for one or multiple OpenFOAM simulation
in a capillary network.
"""

import argparse
import matplotlib.pyplot as plt

from HbO2.plot import styles
from HbO2.plot.hemoglobingraph import HemoglobinOnSegmentsPlotter
from HbO2.postprocess.factory.case import make_post_processor
from plot.figureoptions import FigureOptions
from plot.utils import add_multi_panel_parser_options, annotate_axis_corner

style_scheme = styles.StyleSchemeCOSH()

parser = argparse.ArgumentParser()
parser.add_argument('--cases', nargs='+', default=['.'],
                    help='Cases for which to plot results')
parser.add_argument('--nProfileMin', '-n', type=int, default=1,
                    help='Minimal number of profiles to plot')
add_multi_panel_parser_options(parser)
figOptions = FigureOptions(parser)
args = parser.parse_args()
cases = args.cases
n_profile_min = args.nProfileMin
figOptions.parseOptions(args)

postprocessors = [make_post_processor(case) for case in cases]
plotters = [HemoglobinOnSegmentsPlotter(pp, n_profile_min=n_profile_min)
            for pp in postprocessors]

styles = [dict() for _ in plotters]
for case_path, style in zip(cases, styles):
    if '20151213' in case_path:
        style.update(style_scheme['MVN1'])
    elif '20150715' in case_path:
        style.update(style_scheme['MVN2'])

if args.multipanel:
    panel_layout = (2, 1) if not args.panelLayout else args.panelLayout
    fig, axs = plt.subplots(*panel_layout, sharex=True)
    plt.sca(axs[0])
    for plotter, style in zip(plotters, styles):
        plotter.plot_std_difference_vs_upstream_std(**style)
    if panel_layout[1] == 1:
        axs[0].set_xlabel('')
    annotate_axis_corner(axs[0], 'A')
    plt.sca(axs[1])
    for plotter, style in zip(plotters, styles):
        plotter.plot_std_slope_vs_upstream_std(**style)
    annotate_axis_corner(axs[1], 'B')
    figOptions.saveFig('plotStdDifferenceAndSlopeVsUpstreamStd')
else:
    plt.clf()
    for plotter, style in zip(plotters, styles):
        plotter.plot_std_difference_vs_upstream_std(**style)
    figOptions.saveFig('plotStdDifferenceVsUpstreamStd')

    plt.clf()
    for plotter, style in zip(plotters, styles):
        plotter.plot_std_slope_vs_upstream_std(**style)
    xlim = (-0.001, plt.gca().get_xlim()[1])
    plt.gca().set_xlim(xlim)
    figOptions.saveFig('plotStdSlopeVsUpstreamStd')

    plt.clf()
    for plotter, style in zip(plotters, styles):
        plotter.plot_std_slope_vs_upstream_mean(**style)
    figOptions.saveFig('plotStdSlopeVsUpstreamMean')

    plt.clf()
    for plotter, style in zip(plotters, styles):
        plotter.plot_std_slope_vs_mean_slope(**style)
    figOptions.saveFig('plotStdSlopeVsMeanSlope')

    plt.clf()
    for plotter, style in zip(plotters, styles):
        plotter.plot_std_slope_vs_length(**style)
    figOptions.saveFig('plotStdSlopeVsLength')
