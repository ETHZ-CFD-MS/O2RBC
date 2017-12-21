"""
Plot the results of a parameter study for PO2 drops
"""

import matplotlib.pyplot as plt
import numpy as np

from HbO2.plot.parameterstudy import ParameterStudyPlotter
from HbO2.plot import labels
from plot.figureoptions import FigureOptions


class PO2DropParamStudyPlotter(ParameterStudyPlotter):

    def __init__(self, post_processor, fig_options):
        """

        Args:
            post_processor (ParameterStudyPostProcessor): parameter study post processor
            fig_options (FigureOptions): figure options
        """
        super(PO2DropParamStudyPlotter, self).__init__(
            post_processor,
            fig_options
        )

    def plot_all(self):
        self.contour_delta_po2()
        self.fig_options.saveFig('contourPlotsDeltaPO2')

    def contour_delta_po2(self):
        levels = np.arange(2, 27, 2)
        vrbc, ld = self.post_processor.param_study.meshGrid()
        delta_p_simul = self._reshape_to_grid(self.post_processor.delta_p_tissue_simul())
        delta_p_fitted = self._reshape_to_grid(self.post_processor.delta_p_tissue_fitted())
        c1 = plt.contour(ld, 1e3*vrbc, delta_p_simul, levels, cmap=plt.cm.jet)
        c2 = plt.contour(ld, 1e3*vrbc, delta_p_fitted, levels, linestyles='dashed', cmap=plt.cm.jet)
        labels.setXLabel('LD', '-')
        labels.setYLabel('vrbc', 'mm/s')
        plt.clabel(c1, inline=1, fontsize=12, fmt='%g')

    def _reshape_to_grid(self, a):
        xx, yy = self.post_processor.param_study.meshGrid()
        return np.reshape(a, xx.shape)
