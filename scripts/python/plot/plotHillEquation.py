#!/usr/bin/env python
"""Plotting functions for the Hill equation."""

import argparse
import matplotlib.pyplot as plt
import numpy as np

from plot.figureoptions import FigureOptions
from HbO2.model.chemistry import BloodChemistry
from HbO2.plot import labels


def plotHillEquation(bloodChem, P, reversed=False):
    S = bloodChem.hillS(P)
    if reversed:
        plt.plot(S, P)
        labels.setXLabel('S', '-')
        labels.setYLabel('P', 'mmHg')
        plt.xlim(xmax=1.0)
    else:
        plt.plot(P, S)
        labels.setXLabel('P', 'mmHg')
        labels.setYLabel('S', '-')
        plt.ylim(ymax=1.0)


def plotdPdS(bloodChem, S):
    dPdS = bloodChem.dPdSHill(S)
    plt.plot(S, dPdS)
    labels.setXLabel('S', '-')
    labels.setYLabel('dP/dS', 'mmHg')
    plt.ylim(ymax=300)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--P50', '-p',   type=float, help='PO2 at half saturation',
                        default=47.9)
    parser.add_argument('-n', type=float, help='Hill exponent', default=2.64)
    figOptions = FigureOptions(parser)
    args  = parser.parse_args()
    P50   = args.P50
    nHill = args.n
    figOptions.parseOptions(args)

    P = np.linspace(0, 100, 200)
    S = np.linspace(0, 0.99, 200)
    bloodChem = BloodChemistry(nHill, P50)

    plotHillEquation(bloodChem, P)
    figOptions.saveFig('hillEquation')
    plt.clf()

    plotHillEquation(bloodChem, P, reversed=True)
    figOptions.saveFig('hillEquationReversed')
    plt.clf()

    plotdPdS(bloodChem, S)
    figOptions.saveFig('hillEquation_dPdS')
    plt.clf()
