"""
Classes for plotting of parameter studies.
"""

from abc import ABCMeta, abstractmethod
import matplotlib.pyplot as plt
import numpy as np

from HbO2.plot import labels


class ParameterStudyPlotter(object):
    """
    Abstract class for plotting a parameter study

    Attributes:
        post_processor: instance of derived class of casePostProcessor
        fig_options: instance of FigureOptions
        param_name_to_locator (dict): dictionary that maps the parameter name to an axis locator
    """
    __metaclass__ = ABCMeta

    def __init__(self, post_processor, fig_options):
        self.post_processor = post_processor
        self.fig_options = fig_options
        self.param_name_to_locator = {}

    @abstractmethod
    def plot_all(self):
        """
        Produces all the plots of the parameter study.
        """
        pass

    def plot_variable_single_parameter(self, values, **kwargs):
        """
        Plot the given values against the (unique) parameter of a parameter study.
        """
        line_style = kwargs.get('style', {'linestyle': '-', 'color': 'k'})
        param_values = np.asarray(
            [p[0] for p in self.post_processor.param_study.paramValues()])
        openfoam_param_name = self.post_processor.param_study.paramPrefixes()[0]
        plotting_param_name = labels.openFOAMVarNameToLatex(openfoam_param_name)
        param_unit = labels.openFOAMVarNameToUnit(openfoam_param_name)
        param_scaling = labels.openFOAMVarScaling(openfoam_param_name)
        plt.plot(param_scaling*param_values, values, **line_style)
        labels.setXLabel(plotting_param_name, param_unit)
        try:
            plt.gca().xaxis.set_major_locator(self.param_name_to_locator[openfoam_param_name])
        except KeyError:
            pass
        plt.xlim((param_scaling*np.min(param_values),
                  param_scaling*np.max(param_values)))
