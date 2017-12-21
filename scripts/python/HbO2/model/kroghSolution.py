#!/usr/bin/env python
"""
Compute Krogh solution for PO2 in tissue, given PO2 of the
capillary.
"""

import numpy as np
from scipy.integrate import ode, quad

from HbO2.model.chemistry import BloodChemistry
from HbO2.model.coefficients import IntravascularResistanceFitter, \
                                    intravascularResistanceFromLDHalf, \
                                    intravascularResistanceAnalytical
from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric


class KroghSolution2DCone(object):
    """Compute the Krogh solution in a two-dimensional cone-shaped tissue domain.

    Calculations based on Luecker et al. (2016) and on further developments.
    """

    # default physiological parameters
    defaultVariables = {
        'NHb': 2.03e1,   # [mol m^-3]
        'VMolO2': 2.544e4, # [mlO2 molO2^-1]
        'RBCVolume': 59e-18, # [m^-3]

        'kappaO2RBC': 9.5e-10,   # [m^2 s^-1]
        'kappaO2Plasma': 2.18e-9,
        'kappaO2Wall': 8.73e-10,
        'kappaO2Tissue': 2.41e-9,

        'alphaRBC': 33.8,        # [mlO2 m^-3 mmHg^-1]
        'alphaPlasma': 28.2,
        'alphaWall': 38.9,
        'alphaTissue': 38.9,

        'P50': 47.9,
        'hillExponent': 2.64,
    }

    def __init__(self, simParams):
        """Instance initialization.

        The default variables, contained in the dictionary defaultVariables, are
        overriden by the corresponding variables in simParams if they are found.

        Args:
            simParams: instance of IOHbO2ParametersAxisymmetric
        """
        self._set_default_variables()
        self.M = simParams['O2ConsumptionRate']
        self.vRBC = simParams['RBCVelocity']
        self.LD = simParams['LDMean']

        # replace default variables by those given in simParams
        for var in self.defaultVariables:
            try:
                setattr(self, var, simParams[var])
            except KeyError:
                pass
        self.chem = BloodChemistry(nHill=self.hillExponent, P50=self.P50)

        # load geometry and compute/read the RBC radius
        self.geometry = simParams.geometry()
        try:
            self.Rrbc = simParams['radiusRBC']
        except (TypeError, KeyError):
            self.Rrbc = 0.75*self.geometry['radiusPlasma']

        # use either PO2RBCInlet or HbInlet for the inlet boundary condition
        try:
            self.PO2Inlet = simParams['PO2RBCInlet']
        except KeyError:
            try:
                self.PO2Inlet = self.chem.hillP(simParams['HbInlet'])
            except KeyError:
                raise KeyError('The key "PO2RBCInlet" or "HbInlet" is required.')

        # fit intravascular resistance coefficient
        ivrFitter = IntravascularResistanceFitter('radiusRBC')
        # ivrFitter = IntravascularResistanceFitter('radiusPlasma')
        self.intravascularResistanceLDHalf = ivrFitter.intravascularResistance(simParams, LD=0.5)
        self.convO2Transport = True

    def _set_default_variables(self):
        # this is not done using setattr so that IDEs understand which attributes the class has
        self.NHb = self.defaultVariables['NHb']
        self.VMolO2 = self.defaultVariables['VMolO2']
        self.RBCVolume = self.defaultVariables['RBCVolume']
        self.kappaO2RBC = self.defaultVariables['kappaO2RBC']
        self.kappaO2Plasma = self.defaultVariables['kappaO2Plasma']
        self.kappaO2Wall = self.defaultVariables['kappaO2Wall']
        self.kappaO2Tissue = self.defaultVariables['kappaO2Tissue']
        self.alphaRBC = self.defaultVariables['alphaRBC']
        self.alphaPlasma = self.defaultVariables['alphaPlasma']
        self.alphaWall = self.defaultVariables['alphaWall']
        self.alphaTissue = self.defaultVariables['alphaTissue']
        self.P50 = self.defaultVariables['P50']
        self.hillExponent = self.defaultVariables['hillExponent']

    def RBCLength(self):
        """Return the RBC length."""
        return self.RBCVolume/(np.pi*self.Rrbc**2)

    def RBCFlux(self):
        return self.LD*self.vRBC/self.RBCLength()

    def RBCO2Capacity(self):
        """Return the oxygen capacity of RBC hemoglobin in mlO2 m^-3."""
        return self.VMolO2*self.NHb

    def RBCConvO2Capacity(self):
        "Return the term Q_maxO2"
        return self.RBCFlux()*self.RBCO2Capacity()*self.RBCVolume

    def totalConvO2Capacity(self, S):
        "Return the term that multiplies dS/dx"
        if self.convO2Transport:
            return self.RBCConvO2Capacity() \
                 + self.alphaEffective()*self.plasmaFlow()*self.chem.dPdSHill(S)
        else:
            return self.RBCConvO2Capacity()

    def plasmaFlow(self):
        Rp = self.geometry['radiusPlasma']
        return self.vRBC*np.pi*Rp**2

    def tubeHematocrit(self):
        Rp = self.geometry['radiusPlasma']
        return self.LD*self.Rrbc**2/Rp**2

    def alphaEffective(self): # effective solubility coefficient in capillary
        Hd = self.tubeHematocrit() # ignore Fahraeus effect
        return Hd*self.alphaRBC + (1. - Hd)*self.alphaPlasma

    def consumptionToX(self, x):
        Rt = self.geometry.tissueRadius(x)
        Ra = self.geometry['radiusTissueLeft']
        Rw = self.geometry['radiusWall']
        meanRtSquared = np.sign(Rt)*1/3.*(Ra**2 + Ra*Rt + Rt**2)
        return x*np.pi*self.M*(meanRtSquared - Rw**2)

    def consumptionPerLengthAtX(self, x):
        Rt = self.geometry.tissueRadius(x)
        Rw = self.geometry['radiusWall']
        RtSquared = np.sign(Rt)*Rt**2
        return np.pi*(RtSquared - Rw**2)*self.M

    def ddxConsumptionPerLengthAtX(self, x):
        Rt = self.geometry.tissueRadius(x)
        Ra = self.geometry['radiusTissueLeft']
        Rv = self.geometry['radiusTissueRight']
        L = self.geometry['domainLength']
        return 2.*np.pi*self.M*Rt*(Rv - Ra)/L

    def saturationAtX(self, x):
        if self.convO2Transport:
            return self.saturationAtXWithO2Convection(x)
        else:
            HbInlet = self.chem.hillS(self.PO2Inlet)
            return HbInlet - self.consumptionToX(x)/self.RBCConvO2Capacity()

    def PO2RBCMeanAtX(self, x):
        return self.chem.hillP(self.saturationAtX(x))

    def odeSourceTerm(self, x):
        sourceTerm = -self.consumptionPerLengthAtX(x)
        if self.convO2Transport:
            Hd = self.tubeHematocrit() # ignore Fahraeus effect
            sourceTerm += self.alphaEffective()*self.plasmaFlow() \
                        * self.intravascResistance()*(1 - Hd) \
                        * self.ddxConsumptionPerLengthAtX(x)
        return sourceTerm

    def odeS(self, x, S):
        Q_O2 = self.totalConvO2Capacity(S)
        return self.odeSourceTerm(x)/Q_O2 if Q_O2 > 0 else -1e10

    def saturationAtXWithO2Convection(self, xValues):
        HbInlet = self.chem.hillS(self.PO2Inlet)
        return self.saturationAtXWithO2ConvectionAndInitialCondition(0, HbInlet, xValues)

    def saturationAtXWithO2ConvectionAndInitialCondition(self, x0, Hb0, xValues):
        r = ode(self.odeS).set_integrator('dopri5', method='bdf', atol=1e-10, rtol=1e-8)
        r.set_initial_value(Hb0, x0).set_f_params()
        try:
            iter(xValues)
        except TypeError:
            return max(r.integrate(xValues)[0], 0.0) if xValues != x0 else Hb0
        else:
            return [max(r.integrate(x)[0], 0.0) if x != x0 else Hb0 for x in xValues]

    def averageSaturation(self):
        L = self.geometry['domainLength']
        x = np.linspace(0, L, 20)
        S = self.saturationAtX(x)
        return np.mean(S)

    def PO2Intravascular(self, xValues):
        S  = self.saturationAtX(xValues)
        Hd = self.tubeHematocrit() # ignore Fahraeus effect
        return self.chem.hillP(S) \
            + (1 - Hd)*self.intravascResistance() \
            * self.consumptionPerLengthAtX(xValues)

    def convectivePO2Drop(self, x):
        return self.PO2Inlet - self.chem.hillP(self.saturationAtX(x))

    def intravascResistance(self):
        return intravascularResistanceFromLDHalf(self.intravascularResistanceLDHalf, self.LD)

    def intravascResistancePO2Drop(self, x):
        return self.intravascResistance()*self.consumptionPerLengthAtX(x)

    def extravascularPO2Drop(self, x, R):
        Rw = self.geometry['radiusWall']
        Rt = self.geometry.tissueRadius(x)
        return -self.M/(4*self.kappaO2Tissue*self.alphaTissue)\
               *(R**2 - Rw**2 - 2*Rt**2*np.log(R/Rw))

    def intraRBCRadialPO2Drop(self, x):
        Rw = self.geometry['radiusWall']
        Rt = self.geometry.tissueRadius(x)
        return self.M/(4*self.LD*self.kappaO2RBC*self.alphaRBC)*(Rt**2 - Rw**2)

    def gradPO2AtRBCMembrane(self, x):
        Rw = self.geometry['radiusWall']
        Rt = self.geometry.tissueRadius(x)
        return -0.5*self.M*(Rt**2 - Rw**2) \
                / (self.LD*self.Rrbc*self.kappaO2RBC*self.alphaRBC)

    # PO2 computation
    def oxygenReleaseRateAnalytical(self, x):
        Rt = self.geometry.tissueRadius(x)
        Rw = self.geometry['radiusWall']
        return self.M/self.LD * (Rt**2 - Rw**2)/self.Rrbc**2

    def PO2RBCAnalytical(self, x, R):
        k = self.oxygenReleaseRateAnalytical(x)
        PO2Centerline = self.PO2Inlet - self.convectivePO2Drop(x) \
                                      + 0.125*k*self.Rrbc**2 \
                                        /(self.kappaO2RBC*self.alphaRBC)
        PO2RBC = PO2Centerline - 0.25*k*R**2/(self.kappaO2RBC*self.alphaRBC)
        return max(PO2RBC, 0.0)

    def PO2PlasmaAnalytical(self, x, R):
        k = self.oxygenReleaseRateAnalytical(x)
        PO2RBCMembrane = self.PO2RBCAnalytical(x, self.Rrbc)
        PO2Plasma = PO2RBCMembrane \
                  - 0.5*k*self.Rrbc**2*(np.log(R/self.Rrbc)/
                                       (self.kappaO2Plasma*self.alphaPlasma))
        return max(PO2Plasma, 0.0)

    def PO2WallAnalytical(self, x, R):
        Rp = self.geometry['radiusPlasma']
        k = self.oxygenReleaseRateAnalytical(x)
        PO2InnerWall = self.PO2PlasmaAnalytical(x, Rp)
        PO2Wall = PO2InnerWall \
                - 0.5*k*self.Rrbc**2*(np.log(R/Rp)
                                     /(self.kappaO2Wall*self.alphaWall))
        return max(PO2Wall, 0.0)

    def PO2TissueAnalytical(self, x, R):
        Rw = self.geometry['radiusWall']
        PO2OuterWall = self.PO2WallAnalytical(x, Rw)
        PO2Tissue = PO2OuterWall - self.extravascularPO2Drop(x, R)
        return max(PO2Tissue, 0.0)

    def PO2Analytical(self, x, R):
        if R <= self.Rrbc:
            return self.PO2RBCAnalytical(x, R)
        elif R <= self.geometry['radiusPlasma']:
            return self.PO2PlasmaAnalytical(x, R)
        elif R <= self.geometry['radiusWall']:
            return self.PO2WallAnalytical(x, R)
        else:
            return self.PO2TissueAnalytical(x, R)

    def PO2Tissue(self, x, R):
        PO2Tissue = self.PO2Inlet - self.convectivePO2Drop(x)    \
                                  - self.intravascResistancePO2Drop(x) \
                                  - self.extravascularPO2Drop(x, R)
        return max(PO2Tissue, 0.0)

    def averagePO2Tissue(self):
        """
        Compute the average value of PO2 in the tissue.
        """
        L = self.geometry['domainLength']
        Ra = self.geometry['radiusTissueLeft']
        Rv = self.geometry['radiusTissueRight']
        Rw = self.geometry['radiusWall']
        slice_integral_f = lambda x: self._sliceIntegralPO2(x)
        integral, error = quad(slice_integral_f, 0, L)
        volume = L*np.pi*(1/3.*(Ra**2 + Ra*Rv + Rv**2) - Rw**2)
        return integral/volume

    def _sliceIntegralPO2(self, x):
        Rw = self.geometry['radiusWall']
        Rt = self.geometry.tissueRadius(x)
        PO2Wall = self.PO2Inlet - self.convectivePO2Drop(x) \
                                - self.intravascResistancePO2Drop(x)
        f = lambda R: self._PO2SliceIntegrand(R, x, PO2Wall)
        integral, error = quad(f, Rw, Rt)
        return integral

    def _PO2SliceIntegrand(self, R, x, PO2Wall):
        return 2*np.pi*R*(PO2Wall - self.extravascularPO2Drop(x, R))

    def PO2TissueArray(self, x, R, LDValues, UValues):
        oldLD = self.LD
        oldVRBC = self.vRBC
        tissuePO2 = np.zeros((len(LDValues), len(UValues)))
        for iLD, LD in enumerate(LDValues):
            for iU, U in enumerate(UValues):
                self.LD = LD
                self.vRBC = U
                tissuePO2[iLD, iU] = self.PO2Tissue(x, R)
        self.LD = oldLD
        self.vRBC = oldVRBC
        return tissuePO2

    def ddxPO2(self, x):
        S = self.saturationAtX(x)
        return self.chem.dPdSHill(S)*self.odeS(x, S)

    ## Functions to evaluate derivatives of PO2 with respect to v_RBC and self.LD

    def convectivePart(self, x):
        "Convective part of the analytical expression for PO2 at the wall"
        P50 = self.chem.P50
        n   = self.chem.nHill
        S0 = self.chem.hillS(self.PO2Inlet)
        J1 = self.consumptionToX(x)*self.RBCLength()/ \
                                    (self.RBCVolume*self.RBCO2Capacity())
        return P50*pow((S0 - J1/(self.vRBC*self.LD))/(1 - S0 + J1/(self.vRBC*self.LD)),1./n)

    def IVResistancePart(self, x):
        """Intravascular resistance part of the analytical expression for PO2
        at the wall"""
        return -self.intravascResistance()*self.consumptionPerLengthAtX(x)

    def ddqRBCConvectivePart(self, x):
        P50 = self.chem.P50
        n   = self.chem.nHill
        S0 = self.chem.hillS(self.PO2Inlet)
        J1 = self.consumptionToX(x)*self.RBCLength()/ \
                                    (self.RBCVolume*self.RBCO2Capacity())
        return P50*J1/n * pow(S0*self.LD*self.vRBC - J1,       1/n - 1) \
                        / pow((1 - S0)*self.LD*self.vRBC + J1, 1/n + 1)

    def ddLDConvectivePart(self, x):
        return self.vRBC*self.ddqRBCConvectivePart(x)

    def ddvRBCConvectivePart(self, x):
        return self.LD*self.ddqRBCConvectivePart(x)

    def ddLDIVResistancePart(self, x):
        return -1./self.LD*self.IVResistancePart(x)

    def ddLDPO2(self, x):
        return self.ddLDConvectivePart(x) + self.ddLDIVResistancePart(x)

    def gradPO2(self, x, LDValues, UValues):
        oldLD = self.LD
        oldVRBC = self.vRBC
        ddLDPO2 = np.zeros((len(LDValues), len(UValues)))
        ddvRBCPO2 = np.zeros((len(LDValues), len(UValues)))
        for iLD, LD in enumerate(LDValues):
            for iU, U in enumerate(UValues):
                self.LD = LD
                self.vRBC = U
                ddLDPO2[iLD, iU]   = self.ddLDConvectivePart(x) + \
                                     self.ddLDIVResistancePart(x)
                ddvRBCPO2[iLD, iU] = self.ddvRBCConvectivePart(x)
        self.LD = oldLD
        self.vRBC = oldVRBC
        return ddLDPO2, ddvRBCPO2

    def normalizedDerivativeRatio(self, x):
        return (self.LD * self.ddLDPO2(x))/(self.vRBC * self.ddvRBCConvectivePart(x))


if __name__ == '__main__':
    simParams = IOHbO2ParametersAxisymmetric('.')
    kroghSol = KroghSolution2DCone(simParams)
    ivrFitter = IntravascularResistanceFitter('radiusRBC')

    vRBC = simParams['RBCVelocity']
    LD = simParams['LDMean']
    L = kroghSol.geometry['domainLength']
    x = 0.5*L
    Rv = kroghSol.geometry['radiusTissueRight']
    R = 17.6e-6

    print "Analytical intravascular resistance for LD={0:f}: {1:g} [mmHg m s/mlO2]" \
            .format(LD, intravascularResistanceAnalytical(simParams))

    print "Interpolated intravascular resistance for LD={0:f}: {1:g} [mmHg m s/mlO2]" \
        .format(LD, ivrFitter.intravascularResistance(simParams))

    print "x = %g" % x
    print "R = %g" % R
    print "LD = %g" % LD
    print "vRBC = %g" % vRBC
    print "PO2 at the inlet:       %f" % kroghSol.PO2Inlet
    print "Hb at the inlet:        %f" % kroghSol.chem.hillS(kroghSol.PO2Inlet)
    print "Saturation at x:        %f" % kroghSol.saturationAtX(x)
    print "Average saturation:     %f" % kroghSol.averageSaturation()
    print "Convective PO2 drop:    %f" % kroghSol.convectivePO2Drop(x)
    print "IV resistance PO2 drop: %f" % kroghSol.intravascResistancePO2Drop(x)
    print "Extravascular PO2 drop: %f" % kroghSol.extravascularPO2Drop(x, R)
    print "Tissue PO2:             %f" % kroghSol.PO2Tissue(x, R)
    print "Average tissue PO2:     %f" % kroghSol.averagePO2Tissue()
    print "Minimal tissue PO2:     %f" % kroghSol.PO2Tissue(L, Rv)
    print "Intra-RBC radial drop:  %f" % kroghSol.intraRBCRadialPO2Drop(x)
    print 'dP/dx at x = 0: {:g}'.format(kroghSol.ddxPO2(0))
    print 'dP/dx at x = L: {:g}'.format(kroghSol.ddxPO2(L))

    PO2Final = 1 + kroghSol.extravascularPO2Drop(L, Rv) + kroghSol.intravascResistancePO2Drop(L)
    HbFinal = kroghSol.chem.hillS(PO2Final)
    hbInlet = kroghSol.saturationAtXWithO2ConvectionAndInitialCondition(L, HbFinal, 0)
    PO2RBCInlet = kroghSol.chem.hillP(hbInlet)
    print hbInlet, PO2RBCInlet
    hbMiddle = kroghSol.chem.hillS(40)
    hbInlet = kroghSol.saturationAtXWithO2ConvectionAndInitialCondition(0.5*L, hbMiddle, 0)
    PO2RBCInlet = kroghSol.chem.hillP(hbInlet)
    print hbInlet, PO2RBCInlet
