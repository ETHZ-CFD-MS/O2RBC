#!/usr/bin/env python
"""
Quantification and fitting of model coefficients for oxygen transport modeling.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import linregress

from HbO2.model.chemistry import BloodChemistry
from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric
from parse.readfile import load_numeric_to_dictionary
from plot.figureoptions import FigureOptions


class IntravascularResistanceFitter(object):
    """Fit the intravascular resistance coefficient based on domain and RBC geometry.

    The intravascular resistance coefficient is read from input files and given
    in units of mmHg m s/mlO2.
    """

    scriptPathEnvVar = 'OF_SCRIPTS'
    radiusRBCFitFilePath = 'input/IVRFittedOnRadiusRBC.txt'
    radiusPlasmaFitFilePath = 'input/IVRFittedOnRadiusPlasma.txt'

    def __init__(self, fitParameter='radiusRBC'):
        self.fitParameter = fitParameter
        if fitParameter not in ['radiusRBC', 'radiusPlasma']:
            raise NameError('Unknown fit parameter {}'.format(fitParameter))
        self._read()

    def intravascularResistance(self, simParams, **kwargs):
        """Return the fitted intravascular resistance in mmHg m s/mlO2.

        The fitting method is chosen based on the attribute fitParameter.

        Args:
            simParams: instance of AxisymmetricSimulationParameters
            LD (optional): if given, overrides simParams['LDMean']
        """
        LD = kwargs.get('LD', simParams['LDMean'])
        fitFunc = getattr(self, 'fit_' + self.fitParameter)
        return fitFunc(simParams, LD)

    def fit_radiusRBC(self, simParams, LD):
        "Fit the IVR using the RBC radius"
        f = interp1d(self.radiusRBCFitValues, self.ivrFittedOnRadiusRBC,
                     fill_value='extrapolate')
        ivrLDHalf = f(simParams[self.fitParameter])
        return intravascularResistanceFromLDHalf(ivrLDHalf, LD)

    def fit_radiusPlasma(self, simParams, LD):
        "Fit the IVR using the plasma radius"
        f = interp1d(self.radiusPlasmaFitValues, self.ivrFittedOnRadiusPlasma,
                     fill_value='extrapolate')
        ivrLDHalf = f(simParams[self.fitParameter])
        return intravascularResistanceFromLDHalf(ivrLDHalf, LD)

    def _read(self):
        """
        Load fitted intravascular coefficients from files.
        """
        scriptPath = os.environ[self.scriptPathEnvVar]
        RBCFilePath = os.path.join(scriptPath, self.radiusRBCFitFilePath)
        plasmaFilePath = os.path.join(scriptPath, self.radiusPlasmaFitFilePath)
        with open(RBCFilePath, 'r') as f:
            data = load_numeric_to_dictionary(f)
            self.radiusRBCFitValues = data['radiusRBC']
            self.ivrFittedOnRadiusRBC = data['IVR']
        with open(plasmaFilePath, 'r') as f:
            data = load_numeric_to_dictionary(f)
            self.radiusPlasmaFitValues = data['radiusPlasma']
            self.ivrFittedOnRadiusPlasma = data['IVR']


def intravascularResistanceFromLDHalf(ivrLDHalf, LD):
    return ivrLDHalf*0.5/LD


def intravascularResistanceAnalytical(simParams):
    "Return the analytical intravascular resistance in mmHg m s/mlO2"
    LD = simParams['LDMean']
    kappaO2RBC = simParams['kappaO2RBC']
    alphaRBC = simParams['alphaRBC']
    kappaO2Plasma = simParams['kappaO2Plasma']
    alphaPlasma = simParams['alphaPlasma']
    kappaO2Wall = simParams['kappaO2Wall']
    alphaWall = simParams['alphaWall']
    Rrbc = simParams['radiusRBC']
    Rp = simParams['radiusPlasma']
    Rw = simParams['radiusWall']
    return 1./(2*np.pi*LD)*(1./(4*kappaO2RBC*alphaRBC)
                            + np.log(Rp/Rrbc)/(kappaO2Plasma*alphaPlasma)
                            + np.log(Rw/Rp  )/(kappaO2Wall  *alphaWall))


class KOSFitter(object):
    """Fits the model coefficient K_OS for COSH."""

    def __init__(self):
        pass

    def fitKOS(self, parameter, KOS):
        x = np.asarray(parameter)
        y = np.asarray(KOS)
        return linregress(x, y)


class DiffusiveInteractionResistanceAnalytical:
    """
    Compute the capillary diffusive interaction resistance analytically based on the linearized model.
    """

    def __init__(self, simParams):
        self.simParams = simParams
        self.ivr_fitter = IntravascularResistanceFitter()

    def k_ci(self):
        """
        Compute the capillary diffusive interaction resistance coefficient.

        Returns:
            (float) diffusive interaction resistance coefficient
        """
        d_t = self.simParams['kappaO2Tissue']
        alpha_t = self.simParams['alphaTissue']
        k_iv = self.ivr_fitter.intravascularResistance(self.simParams)
        L = self.simParams.geometry()['domainLength']
        r_t_mean = self.simParams.geometry().tissueRadius(x=L/2)
        r_w = self.simParams['radiusWall']
        return k_iv + 1/(2*np.pi*d_t*alpha_t)*(np.log(r_t_mean/r_w) - 0.5)

    def decay_time(self, ld, hb):
        """
        Compute the decay time for capillary diffusive interaction

        Args:
            ld (float): linear density
            hb (float): mean hemoglobin saturation

        Returns:
            (float) decay time
        """
        d_t = self.simParams['kappaO2Tissue']
        alpha_t = self.simParams['alphaTissue']
        r_rbc = self.simParams['radiusRBC']
        r_w = self.simParams['radiusWall']
        L = self.simParams.geometry()['domainLength']
        r_t_mean = self.simParams.geometry().tissueRadius(x=L/2)
        n_hb = self.simParams['NHb']
        v_mol_o2 = self.simParams['VMolO2']
        ivr_ld_half = self.ivr_fitter.intravascularResistance(self.simParams, LD=0.5)
        chem = BloodChemistry(nHill=self.simParams['hillExponent'], P50=self.simParams['P50'])

        slice_o2 = np.pi*r_rbc**2*n_hb*v_mol_o2   # mlO2 m^-1
        ld_term = ld/(2*np.pi*d_t*alpha_t)*(np.log(r_t_mean/r_w) - 0.5)  # mmHg m s mlO2^-1
        dp_ds = chem.dPdSHill(hb) # mmHg
        return slice_o2*(0.5*ivr_ld_half + ld_term)/dp_ds


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    figOptions = FigureOptions(parser)
    args = parser.parse_args()
    figOptions.parseOptions(args)

    ivrFitter = IntravascularResistanceFitter()
    plt.plot(1e6*ivrFitter.radiusRBCFitValues, 1e-6*ivrFitter.ivrFittedOnRadiusRBC)
    plt.xlabel(r'$r_c\, [\mu m]$')
    plt.ylabel(r'$K_{IV}$')
    figOptions.saveFig('plot_KIV_vs_rRBC')
    plt.clf()

    plt.plot(1e6*ivrFitter.radiusPlasmaFitValues, 1e-6*ivrFitter.ivrFittedOnRadiusPlasma)
    plt.xlabel(r'$r_p\, [\mu m]$')
    plt.ylabel(r'$K_{IV}$')
    figOptions.saveFig('plot_KIV_vs_rPlasma')

    sim_params = IOHbO2ParametersAxisymmetric('.')
    diff_interaction_coeffs = DiffusiveInteractionResistanceAnalytical(sim_params)
    print diff_interaction_coeffs.decay_time(0.6, 0.65)