#!/usr/bin/env python
#
# Plots results from plotPO2DerivativeRatioAnalytical for several values of (LD, v_RBC).
#

import argparse
import matplotlib.pyplot as plt

from HbO2.plot.analytical import PO2DerivativeRatioPlotter
from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric

from plot.figureoptions import FigureOptions

if __name__ == '__main__':

    # parse arguments
    parser = argparse.ArgumentParser()
    figOptions = FigureOptions(parser)

    args  = parser.parse_args()
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    simParams = IOHbO2ParametersAxisymmetric('.')

    plotter = PO2DerivativeRatioPlotter(simParams)
    LDValues = [0.35, 0.2, 0.5]
    UValues  = [1.2e-3, 2e-3, 0.6e-3]

    fig, axs = plt.subplots(3, 1)
    # fig, axs = plt.subplots(1,3)

    for (ax, LD, U) in zip(axs, LDValues, UValues):
        plotter.kroghSol.LD = LD
        plotter.kroghSol.vRBC = U
        plt.sca(ax)
        plotter.plotPO2DerivativeRatio()
        xlabel      = ax.get_xlabel()
        xticklabels = ax.get_xticks().tolist()
        xticklabels = ['%i' % label for label in xticklabels]
        ax.set_xlabel('')
        ax.set_xticklabels([''])
        # ylabel      = ax.get_ylabel()
        # yticklabels = ax.get_yticks().tolist()
        # yticklabels = ['%i' % label for label in yticklabels]
        # ax.set_ylabel('')
        # ax.set_yticklabels([''])

    axs[-1].set_xlabel(xlabel)
    axs[-1].set_xticklabels(xticklabels)
    # axs[0].set_ylabel(ylabel)
    # axs[0].set_yticklabels(yticklabels)

    figOptions.adjustAxes()

    # put A, B, C at the top right corner of each plot
    axs[0].annotate(r'$\mathbf{A}$', xy=(0.90, 0.91), xycoords='axes fraction')
    axs[1].annotate(r'$\mathbf{B}$', xy=(0.90, 0.91), xycoords='axes fraction')
    axs[2].annotate(r'$\mathbf{C}$', xy=(0.90, 0.91), xycoords='axes fraction')

    plotName = 'plotMultiPO2DerivativeRatio'
    # plotName = 'plotMultiPO2DerivativeRatioHorizontal'
    figOptions.saveFig(plotName)
