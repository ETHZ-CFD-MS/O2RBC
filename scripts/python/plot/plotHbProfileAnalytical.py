#!/usr/bin/env python
"""
Plot Hb profiles in a capillary using the analytical model
based on Krogh's assumptions
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from HbO2.model.kroghSolution import KroghSolution2DCone
from HbO2.parse import LDRBCVelocityOptions
from HbO2.plot import labels
from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric
from plot.figureoptions import FigureOptions


def plotHbProfile(simParams, **kwargs):
    kroghSol = KroghSolution2DCone(simParams)
    kroghSol.intravascularResistanceLDHalf = kwargs.get('K0', kroghSol.intravascularResistanceLDHalf)
    style = kwargs.get('style', {'color': 'k'})

    L = kroghSol.geometry['domainLength']
    nx = 100
    xValues = np.linspace(0, L, nx)

    # kroghSol.convO2Transport = False
    # HbSimplified = kroghSol.saturationAtX(xValues, LD, vRBC)
    # print "S at x = %g: %g (simplified)" % (x, kroghSol.saturationAtX(x, LD, vRBC))

    kroghSol.convO2Transport = True
    Hb = kroghSol.saturationAtXWithO2Convection(xValues)
    x = L
    print "S at x = %g: %g" % (x, kroghSol.saturationAtX(x))

    # plt.plot(1e6*xValues, HbSimplified, color='r', linestyle='-')
    plt.plot(1e6*xValues, Hb, linewidth=1, **style)
    plt.xlim(0, 1e6*simParams['domainLength'])
    labels.setXLabel('x', 'um')
    labels.setYLabel('S')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    LDRBCVelocityOptions.addOptions(parser)
    figOptions = FigureOptions(parser)

    args = parser.parse_args()
    figOptions.parseOptions(args)

    simParams = IOHbO2ParametersAxisymmetric('.')

    plotHbProfile(simParams)
    plotName = 'HbProfileAnalytical'

    figOptions.saveFig(plotName)




