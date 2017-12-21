"""
Postprocessing for simulations that quantify the relation between plasma oxygenation
and RBC saturation/oxygenation
"""

from collections import OrderedDict
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize.minpack import curve_fit

from HbO2.postprocess.po2drop import PO2DropPostProcessor
from HbO2.postprocess.parameterstudy import ParameterStudyPostProcessorDecorator

from postprocessing.sampledsets import SampledSet
from utilities.decorators import lazy_function


class PlasmaOxygenationPostProcessor(PO2DropPostProcessor):
    """
    Postprocessing for simulations that quantify the relation between plasma oxygenation
    and RBC saturation/oxygenation
    """

    result_output_dict = OrderedDict()
    result_output_dict['Hb'] = ('hb', (), '{:7.5g}')
    result_output_dict['PO2Plasma'] = ('po2_plasma', (), '{:7.5g}')
    result_output_dict['PeqHb'] = ('p_eq_hb', (), '{:7.5g}')
    result_output_dict['PO2RBC'] = ('po2_rbc', (), '{:7.5g}')
    result_output_dict['PWall'] = ('p_wall', (), '{:7.5g}')
    result_output_dict['DeltaPeq'] = ('delta_p_eq', (), '{:7.5g}')
    result_output_dict['DeltaPO2RBC'] = ('delta_pO2_rbc', (), '{:7.5g}')

    output_file_name = 'plasmaOxygenationResults.txt'

    def __init__(self, decorated, **kwargs):
        super(PlasmaOxygenationPostProcessor, self).__init__(decorated, **kwargs)
        self.po2_plasma_sampled_dir_name = kwargs['PO2PlasmaSampledSetDir']
        self.po2_plasma_sampled_field = kwargs['PO2PlasmaSampledField']
        self.sampled_set_radii = kwargs['sampledSetRadii']
        self.sampled_set_names = ['r_{:g}_um'.format(r).replace('.', '_')
                                  for r in self.sampled_set_radii]

    @lazy_function
    def po2_plasma(self):
        """
        Compute the time-averaged plasma PO2 at x. An average in cylindrical coordinates is computed
        from the sampled sets available in the plasma.

        Returns:
            float, plasma PO2
        """
        sampled_set = SampledSet(self.case_path, self.po2_plasma_sampled_dir_name)
        x_values = []
        po2_values = []
        for set_name in self.sampled_set_names:
            values = sampled_set.last_time_values(set_name, self.po2_plasma_sampled_field)
            x_values = values[:,0]
            po2_values.append(values[:,1])
        po2_values = sum([r*p_vals for r, p_vals in zip(self.sampled_set_radii, po2_values)])\
                     /sum(self.sampled_set_radii)
        f = interp1d(x_values, po2_values)
        return float(f(self.x))

    def delta_p_eq(self):
        return self.p_eq_hb() - self.po2_plasma()

    def delta_pO2_rbc(self):
        return self.po2_rbc() - self.po2_plasma()

    def oxygen_extraction(self):
        return self._krogh_sol.consumptionPerLengthAtX(self.x)


class LDjtPlasmaOxygenationParameterStudyPostProcessor(ParameterStudyPostProcessorDecorator):

    def __init__(self, decorated, **kwargs):
        super(LDjtPlasmaOxygenationParameterStudyPostProcessor, self).__init__(decorated)
        # self.fit_func = self.fit_func_norm_spacing
        # self.fit_func = self.fit_func_ld_inverse
        # self.fit_func = self.fit_func_weighted_average
        # self.fit_func = self.fit_func_weighted_average2
        self.fit_func = self.fit_func_semi_empirical

    def run(self):
        super(LDjtPlasmaOxygenationParameterStudyPostProcessor, self).run()
        print 'Fit for plasma oxygenation as a function of LD and jt:'
        print 'Fitted coefficient:', self.k_fit()
        print 'Std. dev. of the fit:', self.std_fit()

    def fit_func_ld_inverse(self, x, K):
        ld = x[0,:]
        j_t = x[1,:]
        return K*j_t/ld

    def fit_func_norm_spacing(self, x, K):
        ld = x[0,:]
        j_t = x[1,:]
        return K*j_t*(1/ld - 1)

    def fit_func_weighted_average(self, x, Kc, Ksl):
        ld = x[0,:]
        j_t = x[1,:]
        return j_t*(Kc*(1/ld - 1) + Ksl/ld)

    def fit_func_weighted_average2(self, x, Kc, Ksl):
        ld = x[0,:]
        j_t = x[1,:]
        simParams = self.param_study.baseCaseSimParams
        geom = simParams.geometry()
        rc = simParams['radiusRBC']
        rp = geom['radiusPlasma']
        return j_t*(Kc*(1/ld - 1) + (rp**2 - rc**2)/rc**2*Ksl/ld)

    def fit_func_semi_empirical(self, x, Kc):
        ld = x[0,:]
        j_t = x[1,:]
        simParams = self.param_study.baseCaseSimParams
        geom = simParams.geometry()
        rc = simParams['radiusRBC']
        rp = geom['radiusPlasma']
        return j_t*Kc*((1/ld - 1) + 0.5*(rp**2 - rc**2)/rc**2/ld)

    def k_fit(self):
        """
        Get the fitted coefficient for the PO2 drop between RBC and plasma.
        """
        print self.fit_delta_p()
        return self.fit_delta_p()[0]

    def std_fit(self):
        """
        Return the standard deviation of the fit
        """
        return np.sqrt(self.fit_delta_p()[1][0])

    def fit_delta_p(self):
        """
        Fit the PO2 difference between RBC and plasma using a fitting function that depends
        on linear density and the oxygen flux out of the capillary (j_t).

        Returns:
            fitted coefficient
        """
        x, delta_p = self._assemble_fit_func_xy()
        return curve_fit(self.fit_func, x, delta_p)

    def delta_p_fitted(self):
        """
        Return the fitted PO2 difference.

        Returns:
            np.ndarray, fitted PO2 difference
        """
        x, delta_p = self._assemble_fit_func_xy()
        return self.fit_func(x, *self.k_fit())

    def delta_p_simul(self):
        """
        Return the simulated PO2 difference.

        Returns:
            np.ndarray, simulated PO2 difference
        """
        return self._assemble_fit_func_xy()[1]

    def oxygen_extraction(self):
        """
        Return the oxygen extraction per length and time for all simulations
        """
        return np.asarray([pp.oxygen_extraction for pp in self.post_processors])

    def m0_from_oxygen_extraction(self, jt):
        """
        Return the metabolic oxygen consumption rate from the oxygen extraction.
        """
        x = self.post_processors[0].x
        geom = self.param_study.baseCaseSimParams.geometry()
        rt = geom.tissueRadius(x)
        rw = geom['radiusWall']
        return jt/(np.pi * (rt**2 - rw**2))


    def _assemble_fit_func_xy(self):
        """
        Assemble the xdata and ydata arrays for curve fitting.

        Returns:
            x (np.ndarray): array with x-values for fitting function
            y (np.ndarray): 1D array with y-values for fitting function

        """
        x = np.zeros((2, self.param_study.nCases()))
        delta_p = np.zeros((self.param_study.nCases(),))
        for i, pp in enumerate(self.post_processors):
            x[0,i] = pp.simParams['LDMean']
            x[1,i] = pp.oxygen_extraction()
            delta_p[i] = pp.delta_p_eq()
        return x, delta_p
