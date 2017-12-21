#!/usr/bin/env python
"""
Classes for modeling of diffusive interactions between capillaries
"""

import os
from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.optimize import brentq

from HbO2.model.coefficients import DiffusiveInteractionResistanceAnalytical, \
                                    intravascularResistanceFromLDHalf
from HbO2.model.kroghSolution import KroghSolution2DCone
from HbO2.setup.geometry import ParallelCapillaries
from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric


class DiffusiveInteractionModel(object):
    """
    Abstract class for modeling of the diffusive interaction between capillaries.
    """
    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def outflux_per_length(self, hb, **kwargs):
        """
        Return a list of oxygen fluxes per length for each capillary.

        Args:
            hb: list with hemoglobin saturation values for each capillary

        Returns:
            a numpy array with the outfluxes per length for each capillary
        """
        pass


class DiffusiveInteractionModelParallel(DiffusiveInteractionModel):
    """
    Abstract class for modeling of the diffusive interaction between parallel
    capillaries.

    Attributes:
        geometries (ParallelCapillaries): geometric data of parallel capillaries
    """

    def __init__(self, parallel_capillaries):
        super(DiffusiveInteractionModelParallel, self).__init__()
        self.geometries = parallel_capillaries

    @abstractmethod
    def outflux_per_length(self, hb, **kwargs):
        super(DiffusiveInteractionModelParallel, self).outflux_per_length(hb)
        if len(hb) != self.geometries.nCapillaries():
            raise ValueError("""The number of elements in hb does not match the
                             number of capillaries.""")


class DiffusiveInteractionModelTwoCapillaries(DiffusiveInteractionModelParallel):
    """
    Abstract class for modeling of the diffusive interaction between two parallel
    capillaries.

    Attributes:
        parallel_capillaries (ParallelCapillaries): geometric data of parallel capillaries
    """

    def __init__(self, parallel_capillaries, casePath='.'):
        super(DiffusiveInteractionModelTwoCapillaries, self).__init__(parallel_capillaries)
        self.simParams = IOHbO2ParametersAxisymmetric(casePath)
        self.kroghSol = KroghSolution2DCone(self.simParams)
        self.k_ci = DiffusiveInteractionResistanceAnalytical(self.simParams).k_ci()

    @abstractmethod
    def outflux_per_length(self, hb, **kwargs):
        super(DiffusiveInteractionModelTwoCapillaries, self).outflux_per_length(hb)
        if len(hb) != 2:
            raise NotImplementedError('Only two capillaries are currently supported')

    def _tissue_radius_from_outflux(self, jt):
        """
        Compute the tissue radius that correspond to a given oxygen outflux.

        Args:
            jt (np.ndarray): oxygen outflux per unit length

        Returns:
            tissue radius
        """
        m_0 = self.simParams['O2ConsumptionRate']
        r_w = self.simParams['radiusWall']
        jt = np.asarray(jt)
        r_w = r_w*np.ones(jt.shape)
        return np.where(jt >= 0, np.sqrt(np.clip(jt, 0, np.inf)/(m_0*np.pi) + r_w**2), r_w)

    def _complementary_oxygen_outflux(self, jt1, x):
        """
        Computes the oxygen outflux jt2 so that jt1 + jt2 compensate the
        oxygen consumption in the current x-slice.

        Args:
            jt1: oxygen outflux out of either capillary
            x: longitudinal position

        Returns:
            oxygen outflux out of the other capillary
        """
        m_0 = self.simParams['O2ConsumptionRate']
        r_w = self.simParams['radiusWall']
        a_tot = self.geometries.totalSliceArea(x)
        return m_0*(a_tot - 2*np.pi*r_w**2) - jt1


class DiffusiveInteractionModelKrogh(DiffusiveInteractionModelTwoCapillaries):
    """
    Models the diffusive interaction between two parallel capillaries using the
    Krogh model.
    """

    model_name = 'krogh'

    def __init__(self, parallel_capillaries, casePath='.'):
        super(DiffusiveInteractionModelKrogh, self).__init__(parallel_capillaries, casePath)

    def outflux_per_length(self, hb, **kwargs):
        super(DiffusiveInteractionModelKrogh, self).outflux_per_length(hb, **kwargs)
        x = kwargs.get('x')
        ld = kwargs.get('ld', self.simParams['LDMean']*np.ones((2,)))
        po2_rbc = [self.kroghSol.chem.hillP(S) for S in hb]
        jt_a = brentq(self._outflux_equation, -1, 1, args=(po2_rbc, x, ld))
        jt_b = self._complementary_oxygen_outflux(jt_a, x)
        return np.array([jt_a, jt_b])

    def _outflux_equation(self, jt_a, po2_rbc, x, ld):
        """
        Equation that should be satisfied by the oxygen flux per unit length
        out of the first capillary.

        Args:
            jt_a: fluxes to the tissue per unit length (list)
            po2_rbc: RBC PO2 values in capillaries (list)
            x: longitudinal position
            ld: linear density (list with length 2)

        Returns:
            Residual
        """
        ivr_ld_half = self.kroghSol.intravascularResistanceLDHalf
        k_iv = [intravascularResistanceFromLDHalf(ivr_ld_half, l) for l in ld]
        jt_b = self._complementary_oxygen_outflux(jt_a, x)
        rt_a = self._tissue_radius_from_outflux(jt_a)
        rt_b = self._tissue_radius_from_outflux(jt_b)
        return (po2_rbc[0] - po2_rbc[1]
                - (k_iv[0]*jt_a - k_iv[1]*jt_b
                + self._extravascular_drop(rt_a)
                - self._extravascular_drop(rt_b)))

    def _extravascular_drop(self, r_t):
        """
        Compute the extravascular PO2 drop as a function of the tissue radius.

        The extravascular drop is calculated using the expression from the Krogh
        solution. If the tissue radius is less than or equal to the wall radius,
        the extravascular drop is set to zero.

        Args:
            r_t: tissue radius

        Returns:
            Matching extravascular drop.
        """
        r_w = self.simParams['radiusWall']
        m_0 = self.simParams['O2ConsumptionRate']
        d_t = self.simParams['kappaO2Tissue']
        alpha_t = self.simParams['alphaTissue']
        # r_t_simplified = self.simParams.geometry().tissueRadius(50e-6)
        if r_t <= r_w:
            return 0
        else:
            # return m_0/(4*d_t*alpha_t)*(2*r_t**2*np.log(r_t_simplified/r_w)
            #                             - r_t**2 + r_w**2)
            return m_0/(4*d_t*alpha_t)*(2*r_t**2*np.log(r_t/r_w)
                                        - r_t**2 + r_w**2)


class DiffusiveInteractionModelKroghSimplified(DiffusiveInteractionModelTwoCapillaries):
    """
    Models the diffusive interaction between parallel capillaries using the
    simplified version of the Krogh model, with an explicit formula for the fluxes
    out of the capillary.
    """

    model_name = 'simple'

    def __init__(self, parallel_capillaries, casePath='.'):
        super(DiffusiveInteractionModelKroghSimplified, self).__init__(parallel_capillaries, casePath)

    def outflux_per_length(self, hb, **kwargs):
        super(DiffusiveInteractionModelKroghSimplified, self).outflux_per_length(hb, **kwargs)
        x = kwargs.get('x')
        po2_rbc = [self.kroghSol.chem.hillP(S) for S in hb]
        a_tot = self.geometries.tissueSliceArea(x)
        m_0 = self.simParams['O2ConsumptionRate']
        jt_a = 0.5*m_0*a_tot \
               + (po2_rbc[0] - po2_rbc[1])/(2 * self.k_ci)
        jt_b = self._complementary_oxygen_outflux(jt_a, x)
        return np.array([jt_a, jt_b])


class DiffusiveInteractionModelKroghLinearized(DiffusiveInteractionModelTwoCapillaries):
    """
    Models the diffusive interaction between parallel capillaries using the
    linearized version of the Krogh model, with an explicit formula for the fluxes
    out of the capillary that is a function of the hemoblogin saturation in both capillaries,
    and of the derivative of the dissociation curve with respect to the hemoglobin saturation.
    """

    model_name = 'linearized'

    def __init__(self, parallel_capillaries, casePath='.'):
        super(DiffusiveInteractionModelKroghLinearized, self).__init__(parallel_capillaries, casePath)

    def outflux_per_length(self, hb, **kwargs):
        super(DiffusiveInteractionModelKroghLinearized, self).outflux_per_length(hb, **kwargs)
        x = kwargs.get('x')
        a_tot = self.geometries.tissueSliceArea(x)
        m_0 = self.simParams['O2ConsumptionRate']
        hb_mean = np.mean(hb)
        dp_ds = self.kroghSol.chem.dPdSHill(hb_mean)
        jt_a = 0.5*m_0*a_tot + 1/(2 * self.k_ci) * (hb[0] - hb[1]) * dp_ds
        jt_b = self._complementary_oxygen_outflux(jt_a, x)
        return np.array([jt_a, jt_b])


class DiffusiveInteractionModelEqualFluxes(DiffusiveInteractionModelTwoCapillaries):
    """
    Models the diffusive interaction between parallel capillaries using the
    the assumption that oxygen outflux out of both capillaries are the same.
    """

    model_name = 'equal_outfluxes'

    def __init__(self, parallel_capillaries, casePath='.'):
        super(DiffusiveInteractionModelEqualFluxes, self).__init__(parallel_capillaries, casePath)

    def outflux_per_length(self, hb, **kwargs):
        super(DiffusiveInteractionModelEqualFluxes, self).outflux_per_length(hb, **kwargs)
        x = kwargs.get('x')
        a_tot = self.geometries.tissueSliceArea(x)
        m_0 = self.simParams['O2ConsumptionRate']
        jt_a = 0.5*m_0*a_tot
        jt_b = jt_a
        return np.array([jt_a, jt_b])


def makeDiffusiveInteractionModelsTwoCapillaries(casePath):
    """
    Construct a list with instances of the concrete models for diffusive interaction
    between two capillaries.

    Args:
        casePath (str): path to case

    Returns:
        Dictionary with instances of subclasses of DiffusiveInteractionModelTwoCapillaries
    """
    diffusive_interaction_models = [DiffusiveInteractionModelKrogh,
                                    DiffusiveInteractionModelKroghSimplified,
                                    DiffusiveInteractionModelKroghLinearized,
                                    DiffusiveInteractionModelEqualFluxes]
    parallel_capillaries = ParallelCapillaries.fromKeyValueFile(
        os.path.join(casePath, 'geometricData'), 2)
    return {cls.model_name: cls(parallel_capillaries, casePath)
            for cls in diffusive_interaction_models}


if __name__ == '__main__':
    parallel_capillaries = ParallelCapillaries.fromJson('geometries.json')
    interaction_model = DiffusiveInteractionModelKrogh(parallel_capillaries)
    fluxes = interaction_model.outflux_per_length([0.7, 0.7], x=100e-6)
    jt_a = fluxes[0]
    jt_b = fluxes[1]
    rt_a = interaction_model._tissue_radius_from_outflux(jt_a)
    rt_b = interaction_model._tissue_radius_from_outflux(jt_b)
    print "Fluxes:       {:g}, {:g}".format(jt_a, jt_b)
    print "Tissue radii: {:g}, {:g}".format(rt_a, rt_b)