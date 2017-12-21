"""
Plot the results of a parameter study with results for OpenFOAM probes.
"""

import matplotlib.pyplot as plt
import numpy as np

from HbO2.plot.parameterstudy import ParameterStudyPlotter
from HbO2.plot import labels
from plot.figureoptions import FigureOptions


class ProbeParamStudyPlotter(ParameterStudyPlotter):

    def __init__(self, post_processor, fig_options):
        """

        Args:
            post_processor (ParameterStudyPostProcessor): parameter study post processor
            fig_options (FigureOptions): figure options
        """
        super(ProbeParamStudyPlotter, self).__init__(
            post_processor,
            fig_options
        )

    def plot_all(self):
        self.contour_time_scale()
        self.fig_options.saveFig('contourPlotTimeScale')

    def contour_time_scale(self):
        levels = np.arange(0, 2, 0.1)
        vrbc, ld = self.post_processor.param_study.meshGrid()
        time_scale = self._reshape_to_grid(self.post_processor.call_post_processor_method(
                                           'exponential_fit_decay_time'))
        c1 = plt.contour(ld, 1e3*vrbc, time_scale, levels, cmap=plt.cm.jet)
        labels.setXLabel('LD', '-')
        labels.setYLabel('vrbc', 'mm/s')
        plt.clabel(c1, inline=1, fontsize=12, fmt='%g')

    def _reshape_to_grid(self, a):
        xx, yy = self.post_processor.param_study.meshGrid()
        return np.reshape(np.asarray(a), xx.shape)
