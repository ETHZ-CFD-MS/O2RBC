"""
Postprocessor for simulations of diffusive interaction between capillaries
and related model.
"""

from collections import OrderedDict
import numpy as np
from HbO2.model.kroghSolution import KroghSolution2DCone

from HbO2.model.coefficients import DiffusiveInteractionResistanceAnalytical
from scipy.optimize.minpack import curve_fit

from HbO2.COSH.integrate import makeDiffusiveInteractionIntegrators
from HbO2.postprocess.case import PostProcessorDecorator


class DiffusiveInteractionPostProcessor(PostProcessorDecorator):
    """
    Postprocessor for simulations of diffusive interaction between capillaries
    and related model.

    Attributes:
        integrators (dict): dictionary with diffusive interaction integrators with model names as keys.
    """

    output_file_name = 'diffusiveInteractionResults.txt'

    result_output_dict = OrderedDict()
    result_output_dict['HbMeanSim'] = ('finalHbMeanFromSimul', (), '{:7.5g}')
    result_output_dict['HbMeanKrogh'] = ('finalHbMeanFromModel', ('krogh',), '{:7.5g}')
    result_output_dict['HbMeanSimple'] = ('finalHbMeanFromModel', ('simple',), '{:7.5g}')
    result_output_dict['absErrorHbMeanKrogh'] = ('absModelErrorInFinalHbMean', ('krogh',), '{:7.5g}')
    result_output_dict['relErrorHbMeanKrogh'] = ('relModelErrorInFinalHbMean', ('krogh',), '{:7.5%}')
    result_output_dict['HbDiffSim'] = ('finalHbDifferenceFromSimul', (), '{:7.5g}')
    result_output_dict['HbDiffKrogh'] = ('finalHbDifferenceFromModel', ('krogh',), '{:7.5g}')
    result_output_dict['HbDiffSimple'] = ('finalHbDifferenceFromModel', ('simple',), '{:7.5g}')
    result_output_dict['HbDiffLinearized'] = ('finalHbDifferenceFromModel', ('linearized',), '{:7.5g}')
    result_output_dict['HbDiffLinearizedODE'] = ('finalHbDifferenceFromLinearizedODE', (), '{:7.5g}')
    result_output_dict['HbDiffEqualOutfluxes'] = ('finalHbDifferenceFromModel', ('equal_outfluxes',), '{:7.5g}')
    result_output_dict['absErrorHbDiffKrogh'] = ('absModelErrorInFinalHbDifference', ('krogh',), '{:7.5g}')
    result_output_dict['absErrorHbDiffSimple'] = ('absModelErrorInFinalHbDifference', ('simple',), '{:7.5g}')
    result_output_dict['absErrorHbDiffLinearized'] = ('absModelErrorInFinalHbDifference', ('linearized',), '{:7.5g}')
    result_output_dict['absErrorHbDiffLinearizedODE'] = ('absLinearizedODEErrorInFinalHbDifference', (), '{:7.5g}')
    result_output_dict['relErrorHbDiffKrogh'] = ('relModelErrorInHbDifferenceDrop', ('krogh',), '{:7.5%}')
    result_output_dict['relErrorHbDiffSimple'] = ('relModelErrorInHbDifferenceDrop', ('simple',), '{:7.5%}')
    result_output_dict['relErrorHbDiffLinearized'] = ('relModelErrorInHbDifferenceDrop', ('linearized',), '{:7.5%}')
    result_output_dict['relErrorHbDiffLinearizedODE'] = ('relLinearizedODEErrorInHbDifferenceDrop', (), '{:7.5%}')
    result_output_dict['simDecayLength'] = ('saturationDifferenceDecayLength', (), '{:7.5g}')
    result_output_dict['theorDecayLength'] = ('theoreticalSaturationDifferenceDecayLength', (), '{:7.5g}')
    result_output_dict['simDecayTime'] = ('saturationDifferenceDecayTime', (), '{:7.5g}')
    result_output_dict['theorDecayTime'] = ('theoreticalSaturationDifferenceDecayTime', (), '{:7.5g}')

    def __init__(self, decorated, **kwargs):
        super(DiffusiveInteractionPostProcessor, self).__init__(decorated)
        self.n_average = kwargs.get('nAverage', 10)
        self.integrators = makeDiffusiveInteractionIntegrators(self.case_path)
        self.geometry = self.simParams.geometry()
        self.xValues = np.linspace(0, self.geometry['domainLength'], 101)

    def first_integrator(self):
        return self.integrators[self.integrators.keys()[0]]

    def first_integrator(self):
        return self.integrators[self.integrators.keys()[0]]

    def initialHbMean(self):
        return np.mean(self.first_integrator().inletHb)

    def initialHbDifference(self):
        return np.abs(np.diff(self.first_integrator().inletHb))[0]

    def initialXCoord(self):
        """
        Return the final x-coordinate for postprocessing.

        The coordinate is chosen so that it is more than half a RBC length
        before the beginning of the computational domain.
        """
        # offset = 5e-6
        # if offset < 0.5*self.simParams['RBCLength']:
        #     warnings.warn('The default offset is too small, replacing it by 0.5*RBCLength.', UserWarning)
        #     offset = 0.5*self.simParams['RBCLength']
        # return offset
        return 0

    def finalXCoord(self):
        """
        Return the final x-coordinate for postprocessing.

        The coordinate is chosen so that it is more than half a RBC length
        before the end of the computational domain.
        """
        return self.simParams['domainLength'] - self.initialXCoord()

    def initialSCoord(self):
        """
        Return the initial s-coordinate for postprocessing.

        The coordinate refers to the same location as x = 0.
        """
        return self.initialXCoord() + self.simParams.sCoordOffset()

    def finalSCoord(self):
        """
        Return the final s-coordinate for postprocessing.

        The coordinate refers to the same location as the one produced by finalXCoord(self).
        """
        return self.finalXCoord() + self.simParams.sCoordOffset()

    def integrator(self, model_name):
        """
        Return the integrator that corresponds to the given model name

        Args:
            model_name (str): model name of the underlying diffusive interaction model
        Returns:
            corresponding integrator
        Raises:
            ValueError if model name does not correspond to any integrator
        """
        try:
            return self.integrators[model_name]
        except KeyError:
            raise ValueError('Invalid model name "{}"'.format(model_name))

    def finalHbFromModel(self, model_name):
        integrator = self.integrator(model_name)
        if self.simParams.cocurrentFlow():
            return integrator.saturationAtX(self.finalXCoord())
        else:
            xValues = np.array([self.initialXCoord(), self.finalXCoord()])
            hb = integrator.saturationAtX(xValues)
            return np.array([hb[0, -1], hb[1, 0]])

    def finalHbMeanFromModel(self, model_name):
        return np.mean(self.finalHbFromModel(model_name))

    def finalHbDifferenceFromModel(self, model_name):
        hb = self.finalHbFromModel(model_name)
        return np.abs(hb[0] - hb[1])

    def finalHbDifferenceFromLinearizedODE(self):
        x_values = np.linspace(self.initialXCoord(), self.finalXCoord(), 100)
        return self.first_integrator().linearizedSaturationDifference(x_values)[-1]

    def finalHbFromSimul(self):
        if self.simParams.cocurrentFlow():
            s = self.finalSCoord()
            hb = [self.rbcDataPostProcessor.fieldAverageOnEdge('Hb_mean', s, i,
                                                               nAverage=self.n_average)
                  for i in range(4)]
        else:
            sCoords = [self.finalSCoord(), self.initialSCoord(),
                       self.initialSCoord(), self.finalSCoord()]
            hb = [self.rbcDataPostProcessor.fieldAverageOnEdge('Hb_mean', s, i,
                                                               nAverage=self.n_average)
                  for i, s in enumerate(sCoords)]
        return np.array(hb)

    def finalHbMeanFromSimul(self):
        hb = self.finalHbFromSimul()
        return np.mean(hb)

    def finalHbDifferenceFromSimul(self):
        hb = self.finalHbFromSimul()
        hb_a = np.mean(hb[np.array([0, 3])])
        hb_b = np.mean(hb[np.array([1, 2])])
        return np.abs(hb_b - hb_a)

    def hbFromSimul(self, x):
        """
        Return the simulated values of hemoglobin saturation on all edges and at the reguired
        x-locations.

        Args:
            x (np.ndarray): x-locations (defined from 0 to the domain length)

        Returns:
            np.ndarray with shape (4, len(x))

        """
        hb = np.zeros((4, len(x)))
        s = x + self.simParams.sCoordOffset()
        for i in range(4):
            hb[i,:] = self.rbcDataPostProcessor.fieldAverageOnEdge('Hb_mean', s, i,
                                                                   nAverage=self.n_average)
        return hb

    def hbMeanFromSimul(self, x):
        hb = self.hbFromSimul(x)
        return np.mean(hb, axis=0)

    def hbDifferenceFromSimul(self, x):
        hb = self.hbFromSimul(x)
        hb_a = 0.5*(hb[0,:] + hb[3,:])
        hb_b = 0.5*(hb[1,:] + hb[2,:])
        return np.abs(hb_b - hb_a)

    def absModelErrorInFinalHbMean(self, model_name):
        simulMean = self.finalHbMeanFromSimul()
        modelMean = self.finalHbMeanFromModel(model_name)
        return modelMean - simulMean

    def relModelErrorInFinalHbMean(self, model_name):
        absError = self.absModelErrorInFinalHbMean(model_name)
        initialMean = self.initialHbMean()
        simulMean = self.finalHbMeanFromSimul()
        return absError/(initialMean - simulMean)

    def absModelErrorInFinalHbDifference(self, model_name):
        simulDiff = self.finalHbDifferenceFromSimul()
        modelDiff = self.finalHbDifferenceFromModel(model_name)
        return modelDiff - simulDiff

    def relModelErrorInHbDifferenceDrop(self, model_name):
        absError = self.absModelErrorInFinalHbDifference(model_name)
        return -absError/self.hbDifferenceDropFromSimul()

    def absLinearizedODEErrorInFinalHbDifference(self):
        simulDiff = self.finalHbDifferenceFromSimul()
        modelDiff = self.finalHbDifferenceFromLinearizedODE()
        return modelDiff - simulDiff

    def relLinearizedODEErrorInHbDifferenceDrop(self):
        absError = self.absLinearizedODEErrorInFinalHbDifference()
        return -absError/self.hbDifferenceDropFromSimul()

    def hbDifferenceDropFromSimul(self):
        return np.abs(self.initialHbDifference() - self.finalHbDifferenceFromSimul())

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

    def fitExponentialToHbDifference(self):
        """
        Fit parameters of an exponential curve to the simulated hemoglobin saturation difference

        Return:
            float, fitted decay length
        """
        diffSim = self.hbDifferenceFromSimul(self.xValues)
        popt, pcov = curve_fit(self.exponentialFitFunc, self.xValues, diffSim, p0=(100e-6, diffSim[0]))
        return popt

    def exponentialFitToHbDifference(self):
        popt = self.fitExponentialToHbDifference()
        return self.exponentialFitFunc(self.xValues, *popt)

    def saturationDifferenceDecayLength(self):
        """
        Compute the decay length of the hemoglobin saturation difference

        Return:
            float, fitted decay length
        """
        return self.fitExponentialToHbDifference()[0]

    def theoreticalSaturationDifferenceDecayLength(self):
        """
        Compute the theoretical decay length of the hemoglobin saturation difference

        Return:
            float, theoretical decay length
        """
        return self.theoreticalSaturationDifferenceDecayTime()*self.simParams['RBCVelocity']

    def saturationDifferenceDecayTime(self):
        """
        Compute the decay time of the hemoglobin saturation difference

        Return:
            float, fitted decay time
        """
        return self.saturationDifferenceDecayLength()/self.simParams['RBCVelocity']

    def theoreticalSaturationDifferenceDecayTime(self):
        """
        Compute the theoretical decay time of the hemoglobin saturation difference

        Return:
            float, theoretical decay time
        """
        ld = self.simParams['LDMean']
        hb_middle = self.hbMeanFromSimul(np.array([0.5*self.geometry['domainLength']]))[0]
        return DiffusiveInteractionResistanceAnalytical(self.simParams).decay_time(ld, hb_middle)
        # hb_mean = np.mean(self.hbMeanFromSimul(self.xValues))
        # return DiffusiveInteractionResistanceAnalytical(self.simParams).decay_time(ld, hb_mean)
        # vrbc = self.simParams['RBCVelocity']
        # kroghSol = KroghSolution2DCone(self.simParams)
        # hb_krogh = kroghSol.saturationAtX(0.5*self.geometry['domainLength'])
        # return DiffusiveInteractionResistanceAnalytical(self.simParams).decay_time(ld, hb_krogh)

