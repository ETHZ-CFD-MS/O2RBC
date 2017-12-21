#!/usr/bin/env python
"""
Plot Hb profiles in a capillary using the analytical model.
One plot in the geometry read from simParams, the other in a
geometry which is a straight cylinder with the same tissue volume.
"""

import argparse

import numpy as np

from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric
from plot.figureoptions import FigureOptions
from plot.plotHbProfileAnalytical import plotHbProfile

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    figOptions = FigureOptions(parser)
    args = parser.parse_args()
    figOptions.parseOptions(args)

    simParams = IOHbO2ParametersAxisymmetric('.')

    taperedStyle = {'style': {'color': 'k', 'linestyle': '-'}}
    straightStyle = {'style': {'color': 'b', 'linestyle': '--'}}
    plotHbProfile(simParams, **taperedStyle)

    Ra = simParams['radiusTissueLeft']
    Rv = simParams['radiusTissueRight']
    equivalentRadius = np.sqrt(1./3.*(Ra**2 + Ra*Rv + Rv**2))
    simParams['radiusTissueLeft'] = equivalentRadius
    simParams['radiusTissueRight'] = equivalentRadius
    plotHbProfile(simParams, **straightStyle)

    figOptions.saveFig('plotHbProfileTaperedAndStraightCylinders')
