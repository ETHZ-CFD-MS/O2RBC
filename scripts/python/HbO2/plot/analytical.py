#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

from HbO2.model.kroghSolution import KroghSolution2DCone


class PO2DerivativeRatioPlotter(object):
    """
    Plot the ratio of normalized PO2 derivatives w.r.t. linear density RBC velocity,
    using the analytical model based on Krogh's assumptions
    """

    hypoxiaThreshold = 2.0  # [mmHg]

    def __init__(self, simParams):
        self.kroghSol = KroghSolution2DCone(simParams)
        self.kroghSol.convO2Transport = False  # for faster computation of PO2Wall
        # self.fig, self.ax = plt.subplots(1)
        self.nx = 50
        self.nM = 50
        # self.fig.canvas.set_window_title('Derivative ratio')

    def plotPO2DerivativeRatio(self):
        ax = plt.gca()
        L = self.kroghSol.geometry['domainLength']
        Rw = self.kroghSol.geometry['radiusWall']
        xValues = np.linspace(L/self.nx, L, self.nx)
        MValues = np.linspace(1, 3000, self.nM)
        xx, MM = np.meshgrid(xValues, MValues)
        derivativeRatio = np.zeros(xx.shape)
        PO2Wall         = np.zeros(xx.shape)
        for j, M in enumerate(MValues):
            self.kroghSol.M = M
            for i, x in enumerate(xValues):
                derivativeRatio[j, i] = self.kroghSol.normalizedDerivativeRatio(x)
                PO2Wall[j,i] = self.kroghSol.PO2Tissue(x, Rw)

        levels = np.arange(1, 6.25, 0.25)

        C1 = ax.contour(1e6*xx, 1e-3*MM, derivativeRatio, levels, zorder=0)
        CLabel = plt.clabel(C1, levels[0:5], inline=True, inline_spacing=40,
                            fontsize=9, fmt='%g', manual=False)
        removeRepeatedLabels(CLabel)

        for l in CLabel:
            l.set_zorder(0)

        ax.contourf(1e6*xx, 1e-3*MM, PO2Wall, [0, self.hypoxiaThreshold],
                         hatches=[None],
                         colors='0.8',
                         zorder=1)

        ax.set_xlim(0, 1e6*max(xValues))
        ax.set_ylim(0, 1e-3*max(MValues))
        ax.set_xlabel(r'$x \; [\muup \mathrm{m}]$')
        ax.set_ylabel(r'$M_0 \; [10^{-3} \mathrm{mlO_2 \, cm^{-3} \, s^{-1}}]$')


def removeRepeatedLabels(labels):
    """Remove repeated texts in a list of text objects."""
    texts = [label.get_text() for label in labels]
    for dup in sorted(list_duplicates(texts)):
        for i in range(1, len(dup)):
            dupIndex = dup[1][i]
            # print "Removing ", labels[dupIndex]
            labels[dupIndex].remove()

def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items()
            if len(locs)>1)
