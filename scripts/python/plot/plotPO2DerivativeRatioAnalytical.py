#!/usr/bin/env python
"""
Plot the ratio of normalized PO2 derivatives w.r.t. linear density RBC velocity,
using the analytical model based on Krogh's assumptions
"""

import argparse

from plot.figureoptions import FigureOptions
from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric
from HbO2.plot.analytical import PO2DerivativeRatioPlotter
from thirdparty.inputexplorer import inputExplorer


def plotPO2DerivativeRatioInteractive(plotter):
    LDinit = 0.3
    vRBCinit = 1.0e-3
    sliders = [{'label': 'LD', 'valmin': 0.1, 'valmax': 0.9, 'valinit': LDinit},
               {'label': 'vRBC', 'valmin': 0.2e-3, 'valmax': 2.4e-3, 'valinit': vRBCinit,
                'valfmt': '%2.4e'}]
    plotter.plotPO2DerivativeRatio(LDinit, vRBCinit)
    inputExplorer(plotter.plotPO2DerivativeRatio, sliders)

if __name__ == '__main__':
    simParams = IOHbO2ParametersAxisymmetric('.')

    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--LD', '-l',   type=float, help='Linear density',
                        default=simParams['LDMean'])
    parser.add_argument('--vRBC', '-v', type=float, help='RBC velocity',
                        default=simParams['RBCVelocity'])
    parser.add_argument('--interactive', '-i', action='store_true', help='Interactive plotting')
    figOptions = FigureOptions(parser)

    args  = parser.parse_args()
    LD          = args.LD
    vRBC        = args.vRBC
    interactive = args.interactive
    figOptions.parseOptions(args)
    simParams['LDMean'] = LD
    simParams['RBCVelocity'] = vRBC

    plotter = PO2DerivativeRatioPlotter(simParams)
    if interactive:
        plotter.nx = 40
        plotter.nM = 40
        plotPO2DerivativeRatioInteractive(plotter)
    else:
        plotter.plotPO2DerivativeRatio()
        figOptions.adjustAxes()
        plotName = 'PO2DerivativeRatio_LD_%g_vRBC_%g' % (LD, vRBC)
        figOptions.saveFig(plotName)

