#!/usr/bin/env python
"""
Plot PO2 profiles in the tissue using the analytical model 
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
from plot.utils import rounded_bounds


def plotPO2Profile(simParams, radius=17.6e-6, **kwargs):
    kroghSol = KroghSolution2DCone(simParams)
    kroghSol.intravascularResistanceLDHalf = kwargs.get('K0', kroghSol.intravascularResistanceLDHalf)
    if 'K0' in kwargs:
        print "plotPO2Profile: Using K0 = %g" % kwargs['K0']
    kroghSol.convO2Transport = kwargs.get('convO2Transport', True)
    style = kwargs.get('style', {'color': 'k'})

    L    = simParams['domainLength']
    nx   = 100
    xValues = np.linspace(0, L, nx)

    PO2 = [kroghSol.PO2Tissue(x, radius) for x in xValues]
    plt.plot(1e6*xValues, PO2, linewidth=1, **style)

    labels.setXLabel('x', 'um')
    labels.setYLabel('PO2', 'mmHg')


def plotPO2ProfileInterstitialSpace(simParams, radius=17.6e-6, **kwargs):
    kroghSol = KroghSolution2DCone(simParams)
    kroghSol.intravascularResistanceLDHalf = kwargs.get('K0', kroghSol.intravascularResistanceLDHalf)
    if 'K0' in kwargs:
        print "plotPO2Profile: Using K0 = %g" % kwargs['K0']
    style = kwargs.get('style', {'color': 'k'})

    L    = simParams['domainLength']
    nx   = 100
    xValues = np.linspace(0, L, nx)

    PO2 = [kroghSol.PO2Tissue(x, radius) for x in xValues]
    plt.plot(1e6*xValues, PO2, linewidth=1, **style)
    plt.ylim(rounded_bounds(PO2, 10))

    interstitialLayerWidth = 0.35e-6
    kroghSol.geometry['radiusWall'] += interstitialLayerWidth
    kroghSol.intravascularResistanceLDHalf *= 1.1
    PO2WithInterstitial = [kroghSol.PO2Tissue(x, radius) for x in xValues]
    plt.plot(1e6*xValues, PO2WithInterstitial, '--')

    labels.setXLabel('x', 'um')
    labels.setYLabel('PO2', 'mmHg')
    print 'Difference: ', np.array(PO2) - np.array(PO2WithInterstitial)


def plotPO2YProfile(simParams, x=150e-6, **kwargs):
    kroghSol = KroghSolution2DCone(simParams)
    kroghSol.intravascularResistanceLDHalf = kwargs.get('K0', kroghSol.intravascularResistanceLDHalf)
    kroghSol.convO2Transport = kwargs.get('convO2Transport', True)
    if 'K0' in kwargs:
        print "plotPO2Profile: Using K0 = %g" % kwargs['K0']
    style = kwargs.get('style', {'color': 'k'})
    PO2AnalyticalStyle = {'color': 'k',
                          'linestyle': '-',
                          'linewidth': 1,
                          'dashes': (1.2,1.5)}
                          # 'dashes': (5,2.5,1,2.5)}
    nr = 300
    rMax = kroghSol.geometry.tissueRadius(x)
    rAnalytical = np.linspace(0, rMax, nr)
    rTissue = np.linspace(kroghSol.geometry['radiusWall'], rMax, nr)

    PO2Analytical = [kroghSol.PO2Analytical(x, r) for r in rAnalytical]
    PO2Tissue = [kroghSol.PO2Tissue(x, r)     for r in rTissue]
    PO2RBCMean = kroghSol.PO2RBCMeanAtX(x)

    plt.plot(1e6*rAnalytical, PO2Analytical, **PO2AnalyticalStyle)
    plt.plot(1e6*rTissue, PO2Tissue, linewidth=1, **style)
    plt.plot(1e6*kroghSol.Rrbc/2, PO2RBCMean, '.', color='k', markersize=4)
    plt.xlim(0, 1e6*rMax)
    plt.ylim(rounded_bounds(PO2Analytical, 10))

    labels.setXLabel('r', 'um')
    labels.setYLabel('PO2', 'mmHg')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    LDRBCVelocityOptions.addOptions(parser)
    figOptions = FigureOptions(parser)
    args = parser.parse_args()
    optionDict = LDRBCVelocityOptions.parseOptions(args)
    figOptions.parseOptions(args)
    
    simParams = IOHbO2ParametersAxisymmetric('.')
    # plotPO2Profile(simParams, optionDict['radius'])
    # plotName = 'PO2ProfileAnalytical_LD_%g_vRBC_%g_R_%g' %\
    #            (simParams['LDMean'], simParams['RBCVelocity'], optionDict['radius'])

    plotPO2ProfileInterstitialSpace(simParams, optionDict['radius'], K0=4.9e6)
    plotName = 'PO2ProfileAnalytical_LD_%g_vRBC_%g_R_%g' % \
               (simParams['LDMean'], simParams['RBCVelocity'], optionDict['radius'])

    # x = 150e-6
    # plotPO2YProfile(simParams, 150e-6)
    # plotName = 'PO2YProfileAnalytical_LD_%g_vRBC_%g_x_%g' % (LD, vRBC, x)

    figOptions.saveFig(plotName)
