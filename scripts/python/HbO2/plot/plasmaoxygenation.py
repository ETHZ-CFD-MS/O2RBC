"""
Plot the results of a parameter study for plasma oxygenation.
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np

from HbO2.plot.parameterstudy import ParameterStudyPlotter
from HbO2.plot import labels
from plot.figureoptions import FigureOptions


class PlasmaOxygenationParamStudyPlotter(ParameterStudyPlotter):

    def __init__(self, post_processor, fig_options):
        """

        Args:
            post_processor (ParameterStudyPostProcessor): parameter study post processor
            fig_options (FigureOptions): figure options
        """
        super(PlasmaOxygenationParamStudyPlotter, self).__init__(
            post_processor,
            fig_options
        )

    def plot_all(self):
        self.contour_delta_po2()
        self.fig_options.saveFig('contourPlotsDeltaPO2')

    def contour_delta_po2(self, manual_labels=False):
        levels = np.concatenate([np.array([2, 4, 6, 8]), np.arange(12, 41, 4)])
        m0, ld = self.post_processor.param_study.meshGrid()
        delta_p_simul = self._reshape_to_grid(self.post_processor.delta_p_simul())
        self.post_processor.fit_func = self.post_processor.fit_func_semi_empirical
        delta_p_fitted_semi_empirical = self._reshape_to_grid(self.post_processor.delta_p_fitted())
        self.post_processor.fit_func = self.post_processor.fit_func_weighted_average
        delta_p_fitted_weighted_average = self._reshape_to_grid(self.post_processor.delta_p_fitted())
        c1 = plt.contour(ld, 1e-3*m0, delta_p_simul, levels, cmap=plt.cm.jet)
        c2 = plt.contour(ld, 1e-3*m0, delta_p_fitted_weighted_average, levels,
                         linestyles='dotted', linewidths=0.8, colors='k')
        c3 = plt.contour(ld, 1e-3*m0, delta_p_fitted_semi_empirical, levels,
                         linestyles='dashed', cmap=plt.cm.jet)
        for c in c2.collections:
            c.set_dashes([(0, (0.8, 2.4))])
        labels.setXLabel('LD', '-')
        labels.setYLabel('M_0', 'M_0')
        if manual_labels:
            fmt={}  # put empty labels in c2 to create space for c1's labels
            for l in c2.levels:
                fmt[l] = '  '
            plt.clabel(c2, inline=True, inline_spacing=26,
                       fontsize=10, fmt=fmt, manual=True)
            plt.clabel(c3, inline=True, inline_spacing=26,
                       fontsize=10, fmt=fmt, manual=True)
            plt.clabel(c1, inline=True, inline_spacing=24,
                       fontsize=10, fmt='%g', manual=True)
        else:
            plt.clabel(c1, inline=1, fontsize=12, fmt='%g')

        # left y-axis
        axleft = plt.gca()
        major_locator = MultipleLocator(0.4)
        major_formatter = FormatStrFormatter('%g')
        minor_locator = MultipleLocator(0.2)
        axleft.yaxis.set_major_locator(major_locator)
        axleft.yaxis.set_major_formatter(major_formatter)
        axleft.yaxis.set_minor_locator(minor_locator)

        # right y-axis
        jt_tick_values = np.linspace(0.5, 3, 6)
        axright = plt.gca().twinx()
        axright.set_ylim(axleft.get_ylim())
        axright.set_yticks(1e3 * self.post_processor.m0_from_oxygen_extraction(1e-12 * jt_tick_values))
        axright.set_yticklabels(['%.1f' % x for x in jt_tick_values])
        labels.setYLabel('j_t', 'j_t')

    def _reshape_to_grid(self, a):
        m0, ld = self.post_processor.param_study.meshGrid()
        return np.reshape(a, m0.shape)
