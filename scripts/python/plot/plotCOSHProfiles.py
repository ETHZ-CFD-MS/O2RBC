#!/usr/bin/env python
"""
Plot longitudinal profiles of Hb with COSH.
"""

import argparse
import matplotlib.pyplot as plt

from HbO2.COSH.integrate import COSHDiscretePDFIntegrator, COSHUniformPDFIntegrator
from HbO2.COSH.plotCOSH import COSHPlotter
from HbO2.postprocess.factory.case import make_post_processor
from HbO2.plot.styles import create_COSH_legend, set_COSH_rc_params
from plot.figureoptions import FigureOptions
from plot.utils import add_multi_panel_parser_options, annotate_axis_corner, MultipleWithMaxNLocator
from utilities.arguments import add_case_argument

set_COSH_rc_params()


def plotHbProfiles(plotter):
    plotter.plotSimMeanHbProfile()
    plotter.plotSimMeanPlusMinusHbProfile()
    if nPlot > 0:
        plotter.plotSimHbProfile(nPlot)
    plotter.plotODEMeanHbProfile()
    plotter.plotODEMeanPlusMinusStdHbProfile()
    plotter.plotNoModelMeanPlusMinusStdHbProfile()
    ax = plt.gca()
    ax.set_xlim(ax.xaxis.get_data_interval())
    ax.yaxis.set_major_locator(MultipleWithMaxNLocator(0.1, 6))


def plotStdHbProfiles(plotter):
    plotter.plotSimStdHbProfile()
    plotter.plotODEStdHbProfile()
    plotter.plotODELinearizedStdHbProfile()
    plotter.plotNoModelStdHbProfile()
    ax = plt.gca()
    create_COSH_legend(ax, bbox_to_anchor=(0.99, 0.82))
    ax.yaxis.set_major_locator(MultipleWithMaxNLocator(0.01, 6))
    ax.set_xlim(ax.xaxis.get_data_interval())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    add_case_argument(parser)
    parser.add_argument('--path', '-p', help='Path of the case to postprocess',
                        default='.')
    parser.add_argument('--nAverage', '-n', type=int, default=10,
                        help='Number of RBC paths used for averaging')
    parser.add_argument('--nPlot', type=int, default=10,
                        help='Number of RBCs used for plotting individual paths')
    add_multi_panel_parser_options(parser)
    figOptions = FigureOptions(parser)
    args = parser.parse_args()
    case_path = args.caseName
    nAverage = args.nAverage
    nPlot = args.nPlot
    figOptions.parseOptions(args)

    postprocessor = make_post_processor(case_path)
    plotter = COSHPlotter(postprocessor)

    if args.multipanel:
        panel_layout = (2, 1) if not args.panelLayout else args.panelLayout
        fig, axs = plt.subplots(*panel_layout, sharex=True)
        plt.sca(axs[0])
        plotHbProfiles(plotter)
        annotate_axis_corner(axs[0], 'A')
        if panel_layout[1] == 1:
            axs[0].set_xlabel('')
        plt.sca(axs[1])
        plotStdHbProfiles(plotter)
        annotate_axis_corner(axs[1], 'B')
        figOptions.saveFig('plotCOSHHbProfilesMultipanel')
    else:
        plotHbProfiles(plotter)
        figOptions.saveFig('plotCOSHHbProfiles')
        plt.figure()
        plotStdHbProfiles(plotter)
        figOptions.saveFig('plotCOSHStdProfiles')
        plt.gca().set_yscale('log')
        figOptions.saveFig('plotCOSHStdProfilesLogScale')

