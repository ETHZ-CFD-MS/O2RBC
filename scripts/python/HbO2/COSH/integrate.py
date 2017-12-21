#!/usr/bin/env python
"""
Implements various integrators for the ODE models for capillary outflow
saturation heterogeneity (COSH)
"""

import copy
import numpy as np
from scipy.integrate import ode
from scipy.interpolate import interp1d

from HbO2.model.diffusiveinteraction import DiffusiveInteractionModelKrogh, \
                                            makeDiffusiveInteractionModelsTwoCapillaries
from HbO2.model.kroghSolution import KroghSolution2DCone
from HbO2.setup.geometry import ParallelCapillaries
from HbO2.setup.rbcPaths import RBCPathGenerator
from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric, \
                                            IOHbO2ParametersStraightCapillaries


class COSHPDFIntegrator(object):
    """Base class for the integration of ordinary differential equations for COSH."""

    def __init__(self, simParams):
        self.simParams = simParams
        self.kroghSol = KroghSolution2DCone(simParams)
        self.KRI = 2*self.kroghSol.intravascResistance()

    @property
    def inletHb(self):
        return self._inletHb

    @inletHb.setter
    def inletHb(self, S):
        self._inletHb = np.asarray(S)

    @property
    def weights(self):
        return self._weights

    @weights.setter
    def weights(self, w):
        self._weights = np.asarray(w)/np.sum(w)

    def weightedMean(self, x):
        """Weighted mean of a 1D array."""
        return np.average(x, weights=self.weights)

    def weightedStd(self, x):
        """Weighted standard deviation of a 1D array."""
        mean = self.weightedMean(x)
        var = np.average((x - mean)**2, weights=self.weights)
        return np.sqrt(var)

    def odeS(self, x, S):
        PEq = self.kroghSol.chem.hillP
        sourceTerm = self.kroghSol.odeSourceTerm(x) \
                     - 1./self.KRI * (PEq(S) - self.weightedMean(PEq(S)))
        return sourceTerm/self.kroghSol.totalConvO2Capacity(S)

    def odeSNoModel(self, x, S):
        return np.array([self.kroghSol.odeS(x, s) for s in S])

    def odeStdSLinearized(self, x, stdS, f_SMean):
        """ODE for the evolution of standard deviation of S.
        
        Args:
            x: longitudinal position
            stdS: function value
            f_SMean: function that computes SMean at x

        Return:
            x-derivative of the standard deviation of hemoglobin saturation
        """
        dPdS = self.kroghSol.chem.dPdSHill
        Q = self.kroghSol.RBCConvO2Capacity()
        return -dPdS(f_SMean(x))/(self.KRI * Q) * stdS

    def saturationAtX(self, xValues):
        """Compute the hemoglobin saturation at the given longitudinal positions.
        If xValues is an array, return a numpy array with shape
        (len(xValues), len(self.inletHb))."""
        r = ode(self.odeS).set_integrator('dopri5', atol=1e-10, rtol=1e-8)
        r.set_initial_value(self.inletHb, xValues[0])
        return integrate(r, xValues, xValues[0], self.inletHb)

    def saturationAtXNoModel(self, xValues):
        """Compute the hemoglobin saturation at the given longitudinal positions.
        If xValues is an array, return a numpy array with shape
        (len(xValues), len(self.inletHb))."""
        r = ode(self.odeSNoModel).set_integrator('dopri5', atol=1e-10, rtol=1e-8)
        r.set_initial_value(self.inletHb, xValues[0])
        return integrate(r, xValues, xValues[0], self.inletHb)

    def linearizedStdSaturation(self, xValues, f_SMean):
        """Compute the linearized standard deviation of hemoglobin saturation.

        Args:
            xValues: longitudinal positions
            f_SMean: function that computes SMean at x

        Return:
            Standard deviation of S at xValues, from the linearized model
        """
        stdS0 = self.weightedStd(self.inletHb)
        r = ode(self.odeStdSLinearized).set_integrator('dopri5',
                                                       atol=1e-10, rtol=1e-8)
        r.set_initial_value(stdS0, xValues[0]).set_f_params(f_SMean)
        return integrate(r, xValues, xValues[0], stdS0)


class COSHDiscretePDFIntegrator(COSHPDFIntegrator):
    """Integrates the ODE model for COSH when the initial Hb distribution
    is a discrete distribution"""

    def __init__(self, simParams):
        """Initialize an instance
        Args:
            simParams: instance of IOHbO2ParametersAxisymmetric
        """
        super(COSHDiscretePDFIntegrator, self).__init__(simParams)
        self.inletHb = map(self.kroghSol.chem.hillS, [simParams['PO2InletLow'], 
                                                      simParams['PO2InletHigh']])
        self.weights = [1.0]*len(self.inletHb)


class COSHUniformPDFIntegrator(COSHPDFIntegrator):
    """Integrates the ODE model for COSH when the initial Hb distribution
    is a discrete distribution"""

    def __init__(self, simParams):
        """Initialize an instance
        Args:
            simParams: instance of IOHbO2ParametersAxisymmetric
        """
        super(COSHUniformPDFIntegrator, self).__init__(simParams)
        self.hb_inlet_min = simParams['HbInletMin']
        self.hb_inlet_max = simParams['HbInletMax']
        self.n_samples = 100
        self.inletHb = np.linspace(self.hb_inlet_min, self.hb_inlet_max, self.n_samples)
        self.weights = [1.0]*len(self.inletHb)


class COSHDiffusiveInteractionIntegrator(object):
    """
    Integrates the hemoglobin saturation equation for two model capillaries with
    diffusive interaction.

    The following use cases are supported:
        - different inlet hemoglobin saturation values
        - different RBC velocities in both capillaries

    Attributes:
        simParams: instance of HbO2ParametersStraightCapillaries
        kroghSol: instance of KroghSolution2DCone
        interactionModel: model for the diffusive interaction (instance
            of DiffusiveInteractionModel)
        inletHb (list): hemoglobin saturation values at the respective capillary inlets.
    """

    def __init__(self, simParams, interactionModel):
        """
        Object initialization.

        Args:
            simParams: instance of IOHbO2ParametersAxisymmetric
            interactionModel: model for the diffusive interaction (instance
                of DiffusiveInteractionModel)
        """
        self.simParams = simParams
        self.kroghSol = KroghSolution2DCone(simParams)
        self.interactionModel = interactionModel
        try:
            self.inletHb = np.array([simParams['HbInlet1'],
                                     simParams['HbInlet2']])
        except KeyError:
            self.inletHb = simParams['HbInlet']*np.ones((2,))
        try:
            self._read_vrbc_from_path_generator()
        except ValueError:
            print "Using the RBC velocity from the SimulationParameters object."
            self.v_rbc = simParams['RBCVelocity']*np.ones((2,))
        try:
            self._read_ld_from_path_generator()
        except ValueError:
            print "Using the linear density from the SimulationParameters object."
            self.ld = simParams['LDMean']*np.ones((2,))

        self.kroghSolCap0 = copy.deepcopy(self.kroghSol)
        self.kroghSolCap0.vRBC = self.v_rbc[0]
        self.kroghSolCap0.LD = self.ld[0]
        self.kroghSolCap1 = copy.deepcopy(self.kroghSol)
        self.kroghSolCap1.vRBC = self.v_rbc[1]
        self.kroghSolCap1.LD = self.ld[1]

    def model_name(self):
        return self.interactionModel.model_name

    def saturationAtX(self, xValues):
        """
        Compute the hemoglobin saturation at the given longitudinal positions.

        Different integration functions are called for cocurrent or countercurrent flow.

        Returns:
            If xValues is a float, returns a float.
            If xValues is an array, return a numpy array with shape
            (len(xValues), len(self.inletHb)).
        """
        if self.simParams.cocurrentFlow():
            return self.saturationAtXCocurrentFlow(xValues)
        else:
            return self.saturationAtXCountercurrentFlow(xValues)

    def saturationAtXCocurrentFlow(self, xValues):
        """
        Compute the hemoglobin saturation at the given longitudinal positions for cocurrent flow.

        Returns:
            If xValues is a float, returns a float.
            If xValues is an array, return a numpy array with shape
            (len(xValues), len(self.inletHb)).
        """
        r = ode(self.odeS_cocurrent_flow).set_integrator\
            ('dopri5', method='bdf', atol=1e-10, rtol=1e-8)
        r.set_initial_value(self.inletHb, 0)
        return integrate(r, xValues, 0, self.inletHb)

    def saturationAtXCountercurrentFlow(self, xValues):
        """
        Compute the hemoglobin saturation at the given longitudinal positions with two
        parallel capillaries that flow in opposite directions.

        We use the convention that the flow direction is positive in the first capillary
        and negative in the second capillary.

        Returns:
            If xValues is a float, returns a float.
            If xValues is an array, return a numpy array with shape
            (len(xValues), len(self.inletHb)).
        """
        r = ode(self.odeS_countercurrent_flow).set_integrator\
            ('dopri5', method='bdf', atol=1e-10, rtol=1e-8)

        # Correction loop for the hemoglobin saturation at x = 0 in the second capillary
        residual = 1e6
        abs_tol = 1e-4
        hb_init_2 = self.predictor_hb_init_countercurrent_flow()
        hb_correction = 0
        while abs(residual) > abs_tol:
            hb_init_2 += hb_correction
            initial_hb = np.array([self.inletHb[0], hb_init_2])
            r.set_initial_value(initial_hb, 0)
            hb_values = integrate(r, xValues, 0, initial_hb)
            #TODO: Make that more robust!
            inlet_hb_2 = hb_values[-1,-1]
            residual = self.simParams['HbInlet2'] - inlet_hb_2
            hb_correction = 0.9*residual
        return hb_values

    def odeS_cocurrent_flow(self, x, S):
        outfluxes = self.interactionModel.outflux_per_length(S, x=x, ld=self.ld)
        return -outfluxes/self.convective_o2_capacities(S)

    def odeS_countercurrent_flow(self, x, S):
        outfluxes = self.interactionModel.outflux_per_length(S, x=x, ld=self.ld)
        sourceTerm = np.array([-1, 1])*outfluxes
        return sourceTerm/self.convective_o2_capacities(S)

    def convective_o2_capacities(self, S):
        return np.array([self.kroghSolCap0.totalConvO2Capacity(S[0]),
                         self.kroghSolCap1.totalConvO2Capacity(S[1])])

    def predictor_hb_init_countercurrent_flow(self):
        """
        Provides a first guess for the hemoglobin saturation at x = 0 in the capillary
        with negative flow direction

        Returns:
            (float)
        """
        L = self.kroghSol.geometry['domainLength']
        hb_drop = self.simParams['HbInlet'] - self.kroghSol.saturationAtX(L)
        return self.inletHb[1] - hb_drop

    def linearizedSaturationDifference(self, xValues):
        """
        Compute the difference in hemoglobin saturation using the linearized model.
        This is based on the assumption that the oxygen convective capacity is the same in
        both model capillaries.

        Args:
            xValues (np.ndarray): values where to integrate.

        Returns:
            np.ndarray with shape of xValues
        """
        inlet_hb_mean = np.mean(self.inletHb)
        inlet_hb_diff = np.abs(np.diff(self.inletHb))
        self.kroghSol.PO2Inlet = self.kroghSol.chem.hillP(inlet_hb_mean)
        hb_mean = self.kroghSol.saturationAtX(xValues)
        f_hb_mean = interp1d(xValues, hb_mean, kind='linear', fill_value=0.5, bounds_error=False)
        r = ode(self.ode_hb_difference_linearized).set_integrator\
            ('dopri5', method='bdf', atol=1e-10, rtol=1e-8)
        r.set_initial_value(inlet_hb_diff, 0).set_f_params(f_hb_mean)
        return integrate(r, xValues, 0, inlet_hb_diff)[:,0]

    def ode_hb_difference_linearized(self, x, delta_hb, f_hb_mean):
        hb_mean = f_hb_mean(x)
        Q = self.kroghSol.totalConvO2Capacity(hb_mean)
        k_ci = self.interactionModel.k_ci
        return -delta_hb/(Q*k_ci)*self.kroghSol.chem.dPdSHill(hb_mean)

    def _read_vrbc_from_path_generator(self):
        rbc_path_generator = RBCPathGenerator(self.simParams.path)
        v_rbc = rbc_path_generator.vRBC
        if v_rbc[0] == v_rbc[3] and v_rbc[1] == v_rbc[2]:
            self.v_rbc = np.array(v_rbc)[:2]
        else:
            raise ValueError('RBC velocities are not the same in diagonally opposed capillaries.')

    def _read_ld_from_path_generator(self):
        rbc_path_generator = RBCPathGenerator(self.simParams.path)
        ld = [distribution['LDMean'] for distribution in rbc_path_generator.LDDistribution]
        if ld[0] == ld[3] and ld[1] == ld[2]:
            self.ld = np.array(ld)[:2]
        else:
            raise ValueError('Linear densities are not the same in diagonally opposed capillaries.')


def integrate(ode, x, x_init, y_init):
    """
    Integrate the given ordinary differential equation on the given values of x with
    the initial condition f(x_init) = y_init.

    Args:
        ode (ode): ode object from scipy.integrate. Initial conditions should be already set
        x (float or iterable): integrand values
        x_init (float): x-value of initial condition
        y_init (float or iterable): y-value of initial condition

    Returns:
            If xValues is a float, returns a float.
            If xValues is an array, return a numpy array with shape (len(x), len(y_init)).

    """
    try:
        iter(x)
    except TypeError:
        return ode.integrate(x) if x!= x_init else y_init
    else:
        return np.asarray([ode.integrate(s) if s != x_init else y_init for s in x])


def makeDiffusiveInteractionIntegrators(casePath):
    """
    Construct a dict with a diffusive interaction integrators for each two-capillary
    diffusive interaction model.

    Args:
        casePath (str): path to case

    Returns:
        Dictionary with integrators. The keys are the diffusive interaction model names.
    """
    simParams = IOHbO2ParametersStraightCapillaries(casePath)
    interaction_models = makeDiffusiveInteractionModelsTwoCapillaries(casePath)
    return {model_name: COSHDiffusiveInteractionIntegrator(simParams, model)
            for model_name, model in interaction_models.iteritems()}


if __name__ == "__main__":
    simParams = IOHbO2ParametersAxisymmetric('.')
    xValues = np.linspace(0, simParams['domainLength'], 40)[1:]
    LD = 0.3
    vRBC = 1.2e-3
    pdfIntegrator = COSHDiscretePDFIntegrator(simParams)
    S = pdfIntegrator.saturationAtX(xValues, LD, vRBC)
    print "Evolution of S for heterogeneous Hb inlet values:"
    print '\n'.join([str(['{:.5f}'.format(x) for x in elem]) for elem in S])

    parallel_capillaries = ParallelCapillaries.fromKeyValueFile('geometricData', 2)
    diffInterModel = DiffusiveInteractionModelKrogh(parallel_capillaries)
    diffInterIntegrator = COSHDiffusiveInteractionIntegrator(simParams,
                                                             diffInterModel)
    S = diffInterIntegrator.saturationAtX(xValues)
    print "Evolution of S in two capillaries with diffusive interaction:"
    print '\n'.join([str(['{:.5f}'.format(x) for x in elem]) for elem in S])
