"""
Computes tissue volumes from data obtained with numerical simulations of oxygen transport.
"""

import numpy as np
from scipy.optimize import brentq

from HbO2.model.kroghSolution import KroghSolution2DCone


class TissueRadiusFitter(object):
    """
    Provides fitting functions for the oxygen extraction from a capillary as a function
    of RBC flow and hemoglobin saturation drop.

    The computations rely on the Krogh cylinder model.

    Attributes:
        sim_params (SimulationParameters): simulation parameters
        krogh_sol (KroghSolution2DCone): solver for Krogh model

    """

    def __init__(self, sim_params):
        """

        Args:
            sim_params (SimulationParameters):
        """
        self.sim_params = sim_params
        self.krogh_sol = KroghSolution2DCone(sim_params)

    def fit_mean_oxygen_extraction_rate(self, hb_a, hb_v):
        """
        Compute the mean oxygen extraction rate that yields the given hemoglobin saturation drop.

        Args:
            hb_a (float): hemoglobin saturation on arterial side
            hb_v (float): hemoglobin saturation on venous side

        Returns:
            float, mean tissue radius
        """
        self.sim_params['HbInlet'] = hb_a
        if abs(self.krogh_sol.RBCFlux() < 1e-12):
            return 0.0
        # save the original method for the oxygen extraction rate to prepare monkey patching
        self.krogh_sol.old_consumptionPerLengthAtX = self.krogh_sol.consumptionPerLengthAtX
        # estimate upper and lower bounds for jt
        max_bound = np.pi*1e4*100e-6**2  # for M0 = 1e4 and r_t = 100um
        min_bound = -max_bound
        # fit jt
        jt = brentq(self._residual_oxygen_extraction_rate_fit, min_bound, max_bound, args=hb_v)
        return jt

    def fit_mean_tissue_radius(self, hb_a, hb_v):
        """
        Compute the mean tissue radius that yields the given hemoglobin saturation drop.

        A tissue radius less than the wall radius can be returned if hb_v is greater than hb_a.

        Args:
            hb_a (float): hemoglobin saturation on arterial side
            hb_v (float): hemoglobin saturation on venous side

        Returns:
            float, mean tissue radius
        """
        jt = self.fit_mean_oxygen_extraction_rate(hb_a, hb_v)
        rt_squared = jt/(self.krogh_sol.M*np.pi) + self.sim_params['radiusWall']**2
        if rt_squared >= 0:
            return np.sqrt(rt_squared)
        else:
            return -np.sqrt(-rt_squared)

    def fit_tissue_volume(self, hb_a, hb_v):
        """
        Compute the tissue volume that yields the given hemoglobin saturation drop.

        Args:
            hb_a (float): hemoglobin saturation on arterial side
            hb_v (float): hemoglobin saturation on venous side

        Returns:
            float, tissue volume
        """
        jt = self.fit_mean_oxygen_extraction_rate(hb_a, hb_v)
        return self.sim_params['domainLength']*jt/self.krogh_sol.M

    def _residual_oxygen_extraction_rate_fit(self, jt, hb_v_goal):
        """
        Compute the deviation between the outflow hemoglobin saturation obtained with
        the given oxygen extraction rate jt and the target value.

        This operation is performed using monkey patching: the method consumptionPerLengthAtX
        of our instance of KroghSolution2DCone is replaced by a method that returns the constant
        value -jt.

        Args:
            jt (float): oxygen extraction rate
            hb_v_goal (float): target outflow hemoglobin saturation

        Returns:
            float, absolute (signed) error
        """
        self.krogh_sol.consumptionPerLengthAtX = lambda x: jt
        hb_v = self.krogh_sol.saturationAtX(self.krogh_sol.geometry['domainLength'])
        return hb_v - hb_v_goal
