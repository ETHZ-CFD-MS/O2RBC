"""
Postprocessing for OpenFOAM probes
"""

from collections import OrderedDict
import numpy as np
from scipy.optimize import curve_fit

from HbO2.postprocess.case import PostProcessorDecorator
from HbO2.postprocess.parameterstudy import ParameterStudyPostProcessorDecorator
from postprocessing.probes import ProbeValue


def exp_fit_func(x, a, b, c):
    return a*np.exp(-x/b) + c


class ProbePostProcessor(PostProcessorDecorator):

    result_output_dict = OrderedDict()
    result_output_dict['expFitTime'] = ('exponential_fit_decay_time', (), '{:7.5g}')
    result_output_dict['expFitAmpl'] = ('exponential_fit_amplitude', (), '{:7.5g}')
    result_output_dict['expFitConst'] = ('exponential_fit_constant', (), '{:7.5g}')
    result_output_dict['maxFitError'] = ('max_fit_error', (), '{:7.5g}')

    output_file_name = 'probeResults.txt'

    def __init__(self, decorated, **kwargs):
        """
        Constructor

        Args:
            decorated (PostProcessDecorator): decorated object
            **probeName (str): name of the probe folder
            **fieldName (str): name of the field
            **probeIndex (int): index of the probe position to fit
            **minTime (float, optional): minimal time for the fit
            **maxTime (float, optional): maximal time for the fit
        """
        super(ProbePostProcessor, self).__init__(decorated)
        self.probe = ProbeValue(self.case_path, kwargs['probeName'], kwargs['fieldName'])
        self.probe_index = kwargs['probeIndex']
        self.fit_min_time = kwargs.get('minTime', -1e300)
        self.fit_max_time = kwargs.get('maxTime', 1e300)

    def exponential_fit(self):
        times = self.probe.times(self.fit_min_time, self.fit_max_time)
        values = self.probe.values(self.probe_index, self.fit_min_time, self.fit_max_time)
        p0 = (values[0] - values[-1], 0.1, values[-1])
        popt, pcov = curve_fit(exp_fit_func, times, values, p0=p0)
        return popt

    def fitted_values(self):
        times = self.probe.times(self.fit_min_time, self.fit_max_time)
        popt = self.exponential_fit()
        return exp_fit_func(times, *popt)

    def exponential_fit_decay_time(self):
        return self.exponential_fit()[1]

    def exponential_fit_constant(self):
        return self.exponential_fit()[2]

    def exponential_fit_amplitude(self):
        return self.exponential_fit()[0]

    def max_fit_error(self):
        values = self.probe.values(self.probe_index, self.fit_min_time, self.fit_max_time)
        return max(abs(self.fitted_values() - values))


class ProbeParameterStudyPostProcessor(ParameterStudyPostProcessorDecorator):

    def __init__(self, decorated, **kwargs):
        super(ProbeParameterStudyPostProcessor, self).__init__(decorated)

    def run(self):
        super(ProbeParameterStudyPostProcessor, self).run()
