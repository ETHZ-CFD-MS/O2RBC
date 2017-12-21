#!/usr/bin/env python
"""
Postprocessing classes for capillary outflow saturation heterogeneity (COSH)
"""

import functools
import os
import sys
import warnings
from collections import OrderedDict

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.optimize.minpack import curve_fit

from HbO2.COSH.integrate import COSHDiscretePDFIntegrator, COSHUniformPDFIntegrator
from HbO2.model.coefficients import IntravascularResistanceFitter, KOSFitter
from HbO2.postprocess import sampledsets
from HbO2.postprocess.case import PostProcessorDecorator
from HbO2.postprocess.parameterstudy import ParameterStudyPostProcessorDecorator
from HbO2.setup.utils import isAxisymmetricCase, isGraphCase
from parse.readfile import load_csv_to_dictionary
from postprocessing.loadSampledSetFiles import loadSampledSets
from utilities.decorators import lazy_property, lazy_function


class COSHPostProcessor(PostProcessorDecorator):
    """
    Postprocessor for simulation of inter-RBC diffusive interaction, with comparison with
    simplified models and parameter fitting.
    """

    sampledDirName = 'yProfiles'
    sampledSetFile = 'midstreamProfile_PO2.xy'

    output_file_name = 'COSHResults.txt'

    oscillationAmplitude = 1.0  # amplitude for the oscillation radius
    relOscillationAmplitude = 0.02  # rel. amplitude for the rel. oscill. rad.

    result_output_dict = OrderedDict()
    result_output_dict['KRI'] = ('KRI', (), '{:7.5g}')
    result_output_dict['KOS'] = ('KOS', (), '{:7.5g}')
    result_output_dict['relErrODE'] = ('relODEModelErrorInStdDrop', (), '{:7.5g}')
    result_output_dict['relErrLin'] = ('relLinearizedModelErrorInStdDrop', (), '{:7.5g}')
    result_output_dict['relErrExpFit'] = ('relExponentialFitErrorInStdDrop', (), '{:7.5g}')
    result_output_dict['decayLength'] = ('stdSaturationDecayLength', (), '{:7.5g}')
    result_output_dict['decayTime'] = ('stdSaturationDecayTime', (), '{:7.5g}')
    result_output_dict['oscillRad'] = ('oscillationRadius', (), '{:7.5g}')
    result_output_dict['relOscillRad'] = ('relativeOscillationRadius', (), '{:7.5g}')
    result_output_dict['intOscillRad'] = ('integralOscillationRadius', (), '{:7.5g}')
    result_output_dict['intOscillDist'] = ('integralOscillationPenetrationDistance', (), '{:7.5g}')
    result_output_dict['relOscAtOscD'] = ('relativeOscillationAtIntegralOscillationRadius', (), '{:7.5g}')
    result_output_dict['normOscillInt'] = ('normalizedTissueOscillationIntegral', (), '{:7.5g}')
    result_output_dict['rModel'] = ('fitRModel', (), '{:7.5g}')

    def __init__(self, decorated, **kwargs):
        """
        Constructors

        Args:
            decorated (PostProcessDecorator): decorated object
            **nAverage (int): number of RBCs used for averaging
        """
        super(COSHPostProcessor, self).__init__(decorated)
        if isAxisymmetricCase(self.case_path):
            self.integrator = COSHDiscretePDFIntegrator(self.simParams)
        elif isGraphCase(self.case_path):
            self.integrator = COSHUniformPDFIntegrator(self.simParams)
        self.nAverage = kwargs['nAverage']
        self.geometry = self.simParams.geometry()
        self.xValues = np.linspace(0, self.geometry['domainLength'], 101)
        try:
            sampledSetDict = loadSampledSets(self.case_path, self.sampledDirName,
                                             self.sampledSetFile,
                                             maxTimeInterval=self.oscillationTimeWindow())
            self.sampledSetStats = sampledsets.SampledSetStatistics(sampledSetDict)
        except IOError:
            warnings.warn("Could not create sampledSetStats.", UserWarning)
            self.sampledSetStats = None

        self.ivrFitter = IntravascularResistanceFitter()

        self.xSampled = 0.5*self.geometry['domainLength']  # position where the radial profiles are sampled
        self.oscillationAmplitude = 1.0  # amplitude for the oscillation radius
        self.relOscillationAmplitude = 0.02  # rel. amplitude for the rel. oscill. rad.

        self.SIntegrated = None  # used to cache integration results
        self.optimBoundsKRI = [1e5, 1e9]

    @lazy_property
    def KRI(self):
        """Get the fitted KRI coefficient."""
        return self.fitKRI()

    @property
    def KOS(self):
        """Get the fitted KOS coefficient."""
        KIV = self.ivrFitter.intravascularResistance(self.simParams)
        return self.KRI - KIV

    def fitKRI(self):
        """Fit the RBC interaction model coefficient to the simulated data"""
        stdSim = self.rbcDataPostProcessor.fieldStd('Hb_mean', self.sValues(),
                                                    nAverage=self.nAverage)
        res = minimize_scalar(self.objectiveFunctionKRI,
                              bounds=self.optimBoundsKRI,
                              args=([stdSim]),
                              method='bounded', options={'disp': True})
        self.checkKRIOptimResult(res)
        return res.x

    def checkKRIOptimResult(self, res):
        """Performs checks on the optimization result for KRI"""
        x = res.x
        relDiffToBounds = [abs(x - c) / c for c in self.optimBoundsKRI]
        if any([d < 1e-3 for d in relDiffToBounds]):
           warnings.warn(
           """The optimization result is very close to the optimization bounds. 
           Check them and rerun the script.""")
        if not res.success:
            warnings.warn(
            "The optimization failed with the following message:\n %s" % res.message)

    def objectiveFunctionKRI(self, KRI, stdSim):
        """Cost function for KRI"""
        self.integrator.KRI = KRI
        # force recomputation of saturationODE
        if hasattr(self, '_lazy_dict_saturationODE'):
            delattr(self, '_lazy_dict_saturationODE')
        self.SIntegrated = self.integrator.saturationAtX(self.xValues)
        stdODE = self.stdSaturationODE()
        return np.sum((stdSim - stdODE)**2)

    @lazy_function
    def saturationODE(self):
        return self.integrator.saturationAtX(self.xValues)

    def saturationNoModel(self):
        return self.integrator.saturationAtXNoModel(self.xValues)

    def meanSaturationODE(self):
        """Return the weighted mean hemoglobin saturation from the ODE model
        at each longitudinal position"""
        return self.meanSaturation(self.saturationODE)

    def meanSaturationNoModel(self):
        """Return the weighted mean hemoglobin saturation with no modeling of RBC interaction
        at each longitudinal position"""
        return self.meanSaturation(self.saturationNoModel)

    def meanSaturation(self, ode):
        """Return the weighted mean hemoglobin saturation from the integration of
        a given ODE at each longitudinal position"""
        weights = self.integrator.weights
        self.SIntegrated = self.integrator.saturationAtX(self.xValues)
        return np.average(ode(), axis=1, weights=weights)

    def stdSaturationODE(self):
        """Return the standard deviation of the hemoglobin saturation
        from the ODE model at each longitudinal position"""
        return self.stdSaturation(self.saturationODE)

    def stdSaturationNoModel(self):
        """Return the standard deviation of the hemoglobin saturation
        without RBC interaction modeling at each longitudinal position"""
        return self.stdSaturation(self.saturationNoModel)

    def stdSaturation(self, ode):
        """Return the standard deviation of the hemoglobin saturation 
        from the integration of a given ode at each longitudinal position"""
        weights = self.integrator.weights
        mean = self.meanSaturation(ode)
        mean = np.array(mean, ndmin=2).transpose()
        diffFromMean = ode() - np.tile(mean, (1, ode().shape[1]))
        var = np.average(diffFromMean**2, axis=1, weights=weights)
        return np.sqrt(var)

    def linearizedStdSaturationODE(self):
        f_SMean = interp1d(self.xValues, self.meanSaturationODE())
        stdS = self.integrator.linearizedStdSaturation(self.xValues, f_SMean)
        return stdS

    def initialStdSaturation(self):
        return self.stdSaturationODE()[0]

    def finalStdSaturationODE(self):
        return self.stdSaturationODE()[-1]

    def finalStdSaturationLinearized(self):
        return self.linearizedStdSaturationODE()[-1]

    def finalStdSaturationExponentialFit(self):
        return self.exponentialFitToStdHb()[-1]

    def finalStdSaturationFromSimul(self):
        return self.rbcDataPostProcessor.fieldStd('Hb_mean', self.sValues()[-1], self.nAverage)

    def stdDropFromSimul(self):
        return self.rbcDataPostProcessor.fieldStd('Hb_mean', self.sValues()[0], self.nAverage) \
               - self.rbcDataPostProcessor.fieldStd('Hb_mean', self.sValues()[-1], self.nAverage)

    def absODEModelErrorInFinalStd(self):
        return self.finalStdSaturationODE() - self.finalStdSaturationFromSimul()

    def absLinearizedModelErrorInFinalStd(self):
        return self.finalStdSaturationLinearized() - self.finalStdSaturationFromSimul()

    def absExponentialFitErrorInFinalStd(self):
        return self.finalStdSaturationExponentialFit() - self.finalStdSaturationFromSimul()

    def relODEModelErrorInStdDrop(self):
        return -self.absODEModelErrorInFinalStd()/self.stdDropFromSimul()

    def relLinearizedModelErrorInStdDrop(self):
        return -self.absLinearizedModelErrorInFinalStd()/self.stdDropFromSimul()

    def relExponentialFitErrorInStdDrop(self):
        return -self.absExponentialFitErrorInFinalStd()/self.stdDropFromSimul()

    @sampledsets.catchAmplitudeError(0, np.nan)
    def oscillationRadius(self):
        """Return the radius where the given amplitude is reached in the tissue."""
        timeWindow = self.oscillationTimeWindow()
        r = self.sampledSetStats.amplitudePosition(self.oscillationAmplitude, timeWindow)
        return r

    @sampledsets.catchAmplitudeError(0, np.nan)
    def relativeOscillationRadius(self):
        """Return the radius where the given relative amplitude is reached 
        in the tissue.
        """
        timeWindow = self.oscillationTimeWindow()
        return self.sampledSetStats.relativeAmplitudePosition(self.relOscillationAmplitude, timeWindow)

    def integralOscillationRadius(self):
        timeWindow = self.oscillationTimeWindow()
        rMin = self.geometry['radiusWall']
        rMax = self.geometry['radiusTissueLeft']
        return self.sampledSetStats.integralOscillationRadius(rMin, rMax, timeWindow)

    def integralOscillationPenetrationDistance(self):
        timeWindow = self.oscillationTimeWindow()
        rMin = self.geometry['radiusWall']
        rMax = self.geometry['radiusTissueLeft']
        return self.sampledSetStats.integralOscillationPenetrationDistance(
                                                            rMin, rMax, timeWindow)

    def relativeOscillationAtIntegralOscillationRadius(self):
        timeWindow = self.oscillationTimeWindow()
        intOscillRad = self.integralOscillationRadius()
        return self.sampledSetStats.interpolatedRelativeAmplitude(intOscillRad, timeWindow)

    def tissueOscillationIntegral(self):
        """Integral of the PO2 oscillation in the tissue."""
        timeWindow = self.oscillationTimeWindow()
        rMin = self.geometry['radiusWall']
        rMax = self.geometry['radiusTissueLeft']
        return self.sampledSetStats.oscillationIntegral(rMin, rMax, timeWindow)

    def normalizedTissueOscillationIntegral(self):
        """Normalized integral of the PO2 oscillation in the tissue."""
        timeWindow = self.oscillationTimeWindow()
        rMin = self.geometry['radiusWall']
        rMax = self.geometry['radiusTissueLeft']
        return self.sampledSetStats.normalizedOscillationIntegral(
                                                        rMin, rMax, timeWindow)

    def oscillationTimeWindow(self):
        """Compute the time window required to compute tissue PO2 oscillations.

        The window is defined as the time required for n RBCs to pass a given point,
        where n is the number of alternating inlet values.
        """
        LD   = self.simParams['LDMean']
        vRBC = self.simParams['RBCVelocity']
        LRBC = self.simParams['RBCLength']
        return len(self.integrator.inletHb) * LRBC / (LD * vRBC)

    @staticmethod
    def exponentialFitFunc(x, x_decay, a):
        """
        Exponential fit function.

        Args:
            x (np.ndarray): function parameter
            x_decay (float): decay length
            a (float): multiplicative factor

        Returns:
            np.ndarray, function result
        """
        return a*np.exp(-x/x_decay)

    def fitExponentialToStdHb(self):
        """
        Fit parameters of an exponential curve to the simulated standard deviation

        Return:
            float, fitted decay length
        """
        stdSim = self.rbcDataPostProcessor.fieldStd('Hb_mean', self.sValues(),
                                                    nAverage=self.nAverage)
        popt, pcov = curve_fit(self.exponentialFitFunc, self.xValues, stdSim, p0=(100e-6, stdSim[0]))
        return popt

    def exponentialFitToStdHb(self):
        popt = self.fitExponentialToStdHb()
        return self.exponentialFitFunc(self.xValues, *popt)

    def stdSaturationDecayLength(self):
        """
        Compute the decay length of the standard deviation of hemoglobin saturation.

        Return:
            float, fitted decay length
        """
        return self.fitExponentialToStdHb()[0]

    def stdSaturationDecayTime(self):
        """
        Compute the decay time of the standard deviation of hemoglobin saturation.

        Return:
            float, fitted decay time
        """
        return self.stdSaturationDecayLength()/self.simParams['RBCVelocity']

    @staticmethod
    def exponentialFitFunc(x, x_decay, a):
        """
        Exponential fit function.

        Args:
            x (np.ndarray): function parameter
            x_decay (float): decay length
            a (float): multiplicative factor

        Returns:
            np.ndarray, function result
        """
        return a*np.exp(-x/x_decay)

    def fitExponentialToStdHb(self):
        """
        Fit parameters of an exponential curve to the simulated standard deviation

        Return:
            float, fitted decay length
        """
        stdSim = self.rbcDataPostProcessor.fieldStd('Hb_mean', self.sValues(),
                                                    nAverage=self.nAverage)
        popt, pcov = curve_fit(self.exponentialFitFunc, self.xValues, stdSim, p0=(100e-6, stdSim[0]))
        return popt

    def exponentialFitToStdHb(self):
        popt = self.fitExponentialToStdHb()
        return self.exponentialFitFunc(self.xValues, *popt)

    def stdSaturationDecayLength(self):
        """
        Compute the decay length of the standard deviation of hemoglobin saturation.

        Return:
            float, fitted decay length
        """
        return self.fitExponentialToStdHb()[0]

    def stdSaturationDecayTime(self):
        """
        Compute the decay time of the standard deviation of hemoglobin saturation.

        Return:
            float, fitted decay time
        """
        return self.stdSaturationDecayLength()/self.simParams['RBCVelocity']

    def fitRModel(self):
        """Fit the model radius so that the analytical KOS fits the simulated KOS."""
        Rw = self.geometry['radiusWall']
        Rt = self.geometry.tissueRadius(self.xSampled)
        f = lambda r: self.analyticalKOS(r, self.xSampled) - self.KOS
        if self.KOS < 0:
            warnings.warn("Negative fitted KOS, cannot fit rModel.")
            rModel = -2
        elif f(Rw)*f(Rt) > 0:
            warnings.warn("rModel cannot be fitted")
            rModel = -1
        else:
            rModel = brentq(f, Rw, Rt)
        return rModel

    def analyticalKOS(self, r, x):
        """Compute the analytical KOS coefficient for given radius and x-position.

        It is based on the extravascular resistance coefficient.
        """
        Rt    = self.geometry.tissueRadius(x)
        Rw    = self.geometry['radiusWall']
        D     = self.simParams['kappaO2Tissue']
        alpha = self.simParams['alphaTissue']
        return 1/(4*np.pi*D*alpha)*((2*Rt**2*np.log(r/Rw) + Rw**2 - r**2)
                                   /(Rt**2 - Rw**2))


class COSHParamStudyPostProcessor(ParameterStudyPostProcessorDecorator):
    """Provides functions to postprocess a parameter study for the study of COSH."""

    output_file_name = 'COSHResults.txt'
    KOSFitFileName = 'KOSFit.txt'

    def __init__(self, decorated, **kwargs):
        super(COSHParamStudyPostProcessor, self).__init__(decorated)
        self.cachePostProcessors = False
        self.KOSFitter = KOSFitter()
        self.initialize_post_processors()
        self.loadPreviousResults()
        for pp in self.post_processors:
            pp.integrator.KRI = pp.KRI
            print "KRI = ", pp.KRI

    def loadPreviousResults(self):
        try:
            with open(os.path.join(self.param_study['path'], self.output_file_name), 'r') as f:
                print "Loading fitted KRI from previous results in {:s}".format(self.output_file_name)
                data = load_csv_to_dictionary(f)
                for i, pp in enumerate(self.post_processors):
                    pp._lazy_KRI = data['KRI'][i]
        except IOError:
            print "No previous results found."

    def run(self):
        self.write_results()
        self.fitKOSWithIntegralOscillationDistance()

    def fitKOSWithIntegralOscillationDistance(self):
        intOscillRad = 1e6*np.asarray([pp.integralOscillationPenetrationDistance()
                                       for pp in self.post_processors])
        KOS = 1e-6*np.asarray([pp.KOS for pp in self.post_processors])
        slope, intercept, rvalue, pvalue, stderr = self.KOSFitter.fitKOS(intOscillRad, KOS)
        with open(self.KOSFitFileName, 'w') as f:
            f.write('Fit of KOS (in mmHg um s/um^3 O2) with the integral oscillation distance (in um):\n')
            f.write('Slope = {}\n'.format(slope))
            f.write('Intercept = {}\n'.format(intercept))
            f.write('r-squared = {}\n'.format(rvalue**2))
            f.write('p-value   = {}\n'.format(pvalue))
        print 'Wrote fitting parameters of KOS to {:s}'.format(self.KOSFitFileName)

    def write_results(self):
        valueHeader = list(self.param_study.paramPrefixes())
        resultHeader = COSHPostProcessor.result_output_dict.keys()
        with open(self.output_file_name, 'w') as f:
            f.write('\t'.join(h for h in valueHeader + resultHeader) + '\n')
            for pp, values in zip(self.post_processors, self.param_study.paramValues()):
                print 'Postprocessing case {:s}...'.format(pp.case_path)
                f.write('\t'.join(['{:g}'.format(v) for v in values]
                                  + [pp.result_string(resultName)
                                     for resultName in resultHeader]))
                f.write('\n')


class COSHLDvRBCParamStudyPostProcessor(COSHParamStudyPostProcessor,
                                        ParameterStudyPostProcessorDecorator):
    """Provides functions to postprocess a LD-U parameter study for COSH."""

    output_file_name = 'COSHResults.txt'

    def __init__(self, decorated, **kwargs):
        super(COSHLDvRBCParamStudyPostProcessor, self).__init__(decorated)
        self.KOSFitFunc = fittingFunctionDistance
        self.KOSFitCoeffs = [121069, 820.917]

    def run(self):
        super(COSHLDvRBCParamStudyPostProcessor, self).run()
        self.fitKOSWithLDAndU()

    def fittedKOS(self):
        """Return the fitted KOS coefficient, without running the optimization."""
        LDValues = [pp.simParams['LDMean'] for pp in self.post_processors]
        UValues  = [pp.simParams['RBCVelocity'] for pp in self.post_processors]
        partialFitFunc = functools.partial(self.KOSFitFunc, self.KOSFitCoeffs)
        return np.asarray(map(partialFitFunc, LDValues, UValues))

    def fitKOSWithLDAndU(self):
        """Fit the model coefficient KOS as a function of LD and vRBC."""
        print "Starting the optimization of the coefficient fits of KOS..."
        res = minimize(self.residualKOSFit, self.KOSFitCoeffs, method='Nelder-Mead')
        print
        print "Optimization completed"
        print res
        print 'Optimization result: ', ', '.join('{:g}'.format(x) for x in res.x)
        self.KOSFitCoeffs = res.x

    def residualKOSFit(self, x):
        """Compute the residual of the fit for KOS.

        Args:
            x: parameters of the fitting function, array-like
        """
        print 'Coefficients: {}\r'.format(x),
        sys.stdout.flush()
        LDValues = [pp.simParams['LDMean'] for pp in self.post_processors]
        UValues  = [pp.simParams['RBCVelocity'] for pp in self.post_processors]
        partialFitFunc = functools.partial(self.KOSFitFunc, x)
        fittedKOS = map(partialFitFunc, LDValues, UValues)
        simulatedKOS = [pp.KOS for pp in self.post_processors]
        # print '\n'.join('{:g} {:g}'.format(a, b) for a, b in zip(fittedKOS, simulatedKOS))
        return sum(np.asarray(fittedKOS) - np.asarray(simulatedKOS))**2

def fittingFunctionInverseFlow(coeffs, LD, U):
    return coeffs[0] + coeffs[1]/(LD*U)

def fittingFunctionDistance(coeffs, LD, U):
    return coeffs[0] + coeffs[1]*(1/LD - 1)/U
