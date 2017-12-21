#!/usr/bin/env python
"""
Plot the sigmoid dissociation curve between oxygen and hemoglobin.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from figureoptions import FigureOptions


def plotDissociationCurve(P50, n):
    x = np.linspace(0, 120, 200)
    y = np.power(x, n)/(np.power(P50, n)+np.power(x, n))

    plt.plot(x, y, 'k-', color='b')
    plt.ylim([0, 1])
    plt.xlabel(r'$\mathrm{PO}_2 \; [\mathrm{mm\,Hg}]$')
    plt.ylabel(r'$\mathrm{S_{O_2}} \; [-]$')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', type=float, help='P_50', default=45)
    parser.add_argument('-n', type=float, help='Hill exponent', default=2.6)
    figOptions = FigureOptions(parser)
    args = parser.parse_args()
    P50 = args.p
    n = args.n
    figOptions.parseOptions(args)

    plotDissociationCurve(P50, n)

    plotName = 'dissociationCurve'
    figOptions.saveFig(plotName)
