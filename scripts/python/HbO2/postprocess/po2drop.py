"""
Postprocessing for the PO2 drop between RBC, wall, tissue etc.
"""

from collections import OrderedDict
import numpy as np
from HbO2.model.kroghSolution import KroghSolution2DCone
from scipy.interpolate import interp1d
from scipy.optimize import leastsq

from HbO2.model.chemistry import BloodChemistry
from HbO2.postprocess.case import PostProcessorDecorator
from HbO2.postprocess.parameterstudy import ParameterStudyPostProcessorDecorator
from postprocessing.sampledsets import SampledSet


class PO2DropPostProcessor(PostProcessorDecorator):
    """
    Postprocessing of the PO2 drop between the RBC, the wall and the tissue.
    """
    result_output_dict = OrderedDict()
    result_output_dict['Hb'] = ('hb', (), '{:7.5g}')
    result_output_dict['PeqHb'] = ('p_eq_hb', (), '{:7.5g}')
    result_output_dict['PWall'] = ('p_wall', (), '{:7.5g}')
    result_output_dict['PTissue'] = ('p_tissue', (), '{:7.5g}')
    result_output_dict['DeltaPWall'] = ('delta_p_wall_simul', (), '{:7.5g}')
    result_output_dict['DeltaPTissue'] = ('delta_p_tissue_simul', (), '{:7.5g}')

    output_file_name = 'PO2DropResults.txt'

    def __init__(self, decorated, **kwargs):
        super(PO2DropPostProcessor, self).__init__(decorated)
        self.sampled_set = SampledSet(self.case_path, kwargs['sampledSetDir'])
        self.wall_set_name = kwargs['wallSetName']
        self.tissue_set_name = kwargs['tissueSetName']
        self.sampled_field = kwargs['sampledField']
        self.tissue_wall_distance = kwargs['tissueWallDistance']
        self.probe_index = kwargs['probeIndex']

        self.r = self.simParams.geometry()['radiusWall'] + kwargs['tissueWallDistance']
        self.x = self.sampled_set.last_time_values(self.tissue_set_name,
                                                   self.sampled_field)[self.probe_index,0]
        self._krogh_sol = KroghSolution2DCone(self.simParams)

    def hb(self):
        return self.rbcDataPostProcessor.fieldAverage('Hb_mean', self.x)

    def po2_rbc(self):
        return self.rbcDataPostProcessor.fieldAverage('PO2_mean', self.x)

    def p_eq_hb(self):
        chem = BloodChemistry(nHill=self.simParams['hillExponent'], P50=self.simParams['P50'])
        return chem.hillP(self.hb())

    def p_wall(self):
        values = self.sampled_set.last_time_values(self.wall_set_name,
                                                   self.sampled_field)
        x = values[:,0]
        po2 = values[:,1]
        f = interp1d(x, po2)
        return float(f(self.x))

    def p_tissue(self):
        return self.sampled_set.last_time_values(self.tissue_set_name,
                                                 self.sampled_field)[self.probe_index,1]

    def delta_p_wall_simul(self):
        return self.p_eq_hb() - self.p_wall()

    def delta_p_tissue_simul(self):
        return self.p_eq_hb() - self.p_tissue()

    def delta_p_wall_fitted(self, ivr_ld_half):
        if not isinstance(ivr_ld_half, np.ndarray):
            ivr_ld_half = np.array([ivr_ld_half])
        self._krogh_sol.intravascularResistanceLDHalf = ivr_ld_half
        return self._krogh_sol.intravascResistancePO2Drop(self.x)[0]

    def delta_p_tissue_fitted(self, ivr_ld_half):
        self._krogh_sol.intravascularResistanceLDHalf = ivr_ld_half
        return self.delta_p_wall_fitted(ivr_ld_half) \
               + self._krogh_sol.extravascularPO2Drop(self.x, self.r)


class PO2DropParameterStudyPostProcessor(ParameterStudyPostProcessorDecorator):

    def __init__(self, decorated, **kwargs):
        super(PO2DropParameterStudyPostProcessor, self).__init__(decorated)

    def run(self):
        super(PO2DropParameterStudyPostProcessor, self).run()
        print 'Fitted IVR with wall PO2:   {:g}'.format(self.k_ivr_fit_with_wall_po2())
        print 'Fitted IVR with tissue PO2:   {:g}'.format(self.k_ivr_fit_with_tissue_po2())

    def k_ivr_fit_with_wall_po2(self):
        """
        Get the fitted IVR coefficient based on wall PO2
        """
        x, cov_x, infodict, mesg, ier = leastsq(self._residual_delta_p_wall, 5e6, xtol=1e-10, ftol=1e-10,
                                                full_output=True)
        if ier not in [1, 2, 3, 4]:
            raise RuntimeError(mesg)
        return x[0]

    def k_ivr_fit_with_tissue_po2(self):
        """
        Get the fitted IVR coefficient based on tissue PO2
        """
        x, cov_x, infodict, mesg, ier = leastsq(self._residual_delta_p_tissue, 5e6, full_output=True)
        if ier not in [1, 2, 3, 4]:
            raise RuntimeError(mesg)
        return x[0]

    def delta_p_tissue_fitted(self):
        ivr_ld_half = self.k_ivr_fit_with_tissue_po2()
        return np.asarray([pp.delta_p_tissue_fitted(ivr_ld_half) for pp in self.post_processors])

    def delta_p_tissue_simul(self):
        return np.asarray([pp.delta_p_tissue_simul() for pp in self.post_processors])

    def _residual_delta_p_wall(self, ivr_ld_half):
        return np.asarray([pp.delta_p_wall_fitted(ivr_ld_half)
                           - pp.delta_p_wall_simul() for pp in self.post_processors])

    def _residual_delta_p_tissue(self, ivr_ld_half):
        return np.asarray([pp.delta_p_tissue_fitted(ivr_ld_half)
                           - pp.delta_p_tissue_simul() for pp in self.post_processors])
