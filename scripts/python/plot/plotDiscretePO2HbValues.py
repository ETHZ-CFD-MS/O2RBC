#!/usr/bin/env python
#
# Plot on the same graph: discrete values of PO2 and hemoblogin
#

import argparse

import matplotlib.pyplot as plt
import numpy as np

from plot.figureoptions import FigureOptions

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    figOptions = FigureOptions(parser)
    args  = parser.parse_args()
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    PO2_values = [16.5, 25.2, 18.6]
    Hb_values  = [0.364, 0.418, 0.363]

    PO2_color = (0.8, 0.0, 0.0)
    Hb_color = (0.0, 0.0, 0.8)

    fig, ax1 = plt.subplots()
    x = range(3)

    ax1.plot(x, PO2_values, 'o', color=PO2_color)
    ax1.set_ylim(10, 30)
    ax1.set_ylabel(r'$\mathrm{PO}_2 \; [\mathrm{mm\,Hg}]$',color=PO2_color)
    for tl in ax1.get_yticklabels():
        tl.set_color(PO2_color)
    
    ax2 = ax1.twinx()
    ax2.plot(x, Hb_values, 's', color=Hb_color)
    ax2.set_ylim(0.35, 0.45)
    ax2.set_yticks(np.arange(0.35, 0.45, 0.025))

    ax2.set_ylabel(r'$\mathrm{S_{O_2} \; [-]$',color=Hb_color)
    for tl in ax2.get_yticklabels():
        tl.set_color(Hb_color)

    # remove xticks and change legend
    ax1.set_xlim(-0.3, 2.3)
    ax1.xaxis.set_ticks(x)
    ax1.tick_params(\
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='on',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='on') # labels along the bottom edge are off

    labels = [item.get_text() for item in ax1.get_xticklabels()]
    labels[0] = r'$\mathrm{Baseline}$' '\n' r'$d = 4\;\mathrm{ \mu m}$' '\n' r'$\mathrm{normal\;CMRO2}$'
    labels[1] = r'$\mathrm{Dilated}$' '\n' r'$d + \mathrm{15\%}$' '\n' r'$\mathrm{normal\;CMRO2}$'
    labels[2] = r'$\mathrm{Dilated}$' '\n' r'$d + \mathrm{15\%}$' '\n' r'$\mathrm{CMRO2 + 15\%}$'

    ax1.set_xticklabels(labels)



    figOptions.adjustAxes()
    figOptions.setGrid()

    plotName = 'dilationValuesPlot'

    figOptions.saveFig(plotName)
