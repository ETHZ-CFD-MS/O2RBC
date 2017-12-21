"""
Classes for postprocessing of parameter studies.
"""

import cPickle as pickle
import csv
import itertools
import os

import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import minimize_scalar

from HbO2.model.kroghSolution import KroghSolution2DCone
from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric
from postprocessing.extractMTC import extractMTCEulerian
from postprocessing.sampledsets import SampledSet


class ParameterStudyPostProcessor(object):
    """
    Class for a parameter study postprocessor
    """

    post_processor_pickle = 'postProcessors.pkl'
    output_file_name = 'parameterStudyResults.txt'

    def __init__(self, param_study, settings_dict):
        self.param_study = param_study
        self.post_processing_settings = settings_dict
        self.post_processors = []
        self.cache_post_processors = False
        self.initialize_post_processors()

    def initialize_post_processors(self):
        if self.cache_post_processors:
            try:
                with open(self.post_processor_pickle, 'rb') as f:
                    print "Reading postprocessors from %s..." % self.post_processor_pickle
                    self.post_processors = pickle.load(f)
                    print "Finished reading postprocessors"
            except IOError:
                print 'Creating postprocessor objects...'
                self.post_processors = map(lambda s: self.make_post_processor(s),
                                           self.param_study.casePaths())
                pickle.dump(self.post_processors, open(self.post_processor_pickle, 'wb'))
                print "Wrote postprocessors to %s" % self.post_processor_pickle
        else:
            self.post_processors = map(lambda s: self.make_post_processor(s),
                                       self.param_study.casePaths())

    def get_post_processors(self):
        """
        Generator for the postprocessors of the individual cases

        Returns:
            Generator
        """
        for case_name in self.param_study.casePaths():
            yield self.make_post_processor(case_name)

    def make_post_processor(self, case_path):
        """
        Abstract factory method for that produces a postprocessor

        Args:
            case_path (str): path to case

        Returns:
           postprocessor object, expected to be a derived class of casePostProcessor
        """
        from HbO2.postprocess.factory.case import make_post_processor
        return make_post_processor(case_path, settings_dict=self.post_processing_settings)

    def run(self):
        pass

    def call_post_processor_method(self, method_name, *args, **kwargs):
        """
        Call a method for each element of the postprocessor sequence and return the result
        in a list.

        Args:
            method_name (str): name of the method to call
            *args: argument sequence to pass to the method
            **args: dictionary of named arguments to pass to the method

        Returns:
            list of result

        Raises:
            AttributeError: if the postprocessors have no method with name method_name
        """
        try:
            return [getattr(pp, method_name)(*args, **kwargs) for pp in self.post_processors]
        except AttributeError:
            print 'No method {:s} for postprocessor instance of type', type(self.post_processors[0])
            raise

    def output_files_and_headers(self):
        """
        Return a dictionary with output file names and keys and list of string for header as value.
        """

        value_header = list(self.param_study.paramPrefixes())
        output_dict = self.post_processors[0].output_files_and_result_names()
        for key in output_dict:
            output_dict[key] = value_header + output_dict[key]
        return output_dict

    def output_files_and_results(self):
        """
        Return a dictionary with output file names and keys and list of results as value.
        """
        output_dict = {}
        for pp, values in zip(self.post_processors, self.param_study.paramValues()):
            pp_dict = pp.output_files_and_result_strings()
            if not output_dict:
                output_dict = {key: [] for key in pp_dict}
            for key in pp_dict:
                param_values = ['{:g}'.format(v) for v in values]
                output_dict[key].append(param_values + pp_dict[key])
        return output_dict


class ParameterStudyPostProcessorDecorator(ParameterStudyPostProcessor):
    """
    Decorator for a parameter study postprocessor.
    """

    def __init__(self, decorated):
        self.decorated = decorated

    def run(self):
        self.decorated.run()

    def __getattr__(self, name):
        return getattr(self.decorated, name)


class ParameterStudyPostProcessorWriter(object):
    """
    Writer for a ParamStudyPostProcessor instance.
    """

    def __init__(self, post_processor):
        self.post_processor = post_processor
        self.separator = '\t'

    def write_results(self):
        files_and_header = self.post_processor.output_files_and_headers()
        files_and_results = self.post_processor.output_files_and_results()
        file_objects = []
        for file_name in files_and_header.keys():
            file_objects.append(open(file_name, 'w'))

        for fo, header, results in zip(file_objects, files_and_header.values(),
                                       files_and_results.values()):
            fo.write(self.separator.join(header))
            fo.write('\n')
            for case_result in results:
                fo.write(self.separator.join(case_result))
                fo.write('\n')
            fo.close()


class LDvRBCParameterStudyPostProcessor(ParameterStudyPostProcessorDecorator):

    PO2MeanOutputFileName = 'PO2MeanParamStudy.txt'
    MTCOutputFileName = 'MTCNuParameterStudy.txt'

    def __init__(self, decorated, **kwargs):
        super(LDvRBCParameterStudyPostProcessor, self).__init__(decorated)
        self.sampledSetDir = kwargs.get('sampledSetDir')
        self.sampledSetName = kwargs.get('sampledSetName')
        self.sampledField = kwargs.get('sampledField')
        self.xFit = kwargs.get('xFit')
        self.rFit = kwargs.get('rFit')
        self.probeSpacing = kwargs.get('probeSpacing')
        convO2Transport = kwargs.get('convO2Transport', True)

        self.kroghSols = [KroghSolution2DCone(IOHbO2ParametersAxisymmetric(p))
                          for p in self.param_study.casePaths()]
        for ks in self.kroghSols:
            ks.convO2Transport = convO2Transport
        simParams = IOHbO2ParametersAxisymmetric(self.param_study.casePaths()[0])
        domainLength = simParams['domainLength']
        self.probeIdx = int(np.round(self.xFit/self.probeSpacing))
        self.probePositions = np.arange(0, domainLength + 1e-6, self.probeSpacing)
        self.probeNames = ['x={:g}um'.format(1e6*x) for x in self.probePositions]

    def run(self):
        # self.write_results()
        self.writePO2Mean()
        # self.writeMTCNu()
        print 'Fitted IVR = {:g}'.format(self.fitIntravascularResistance())

    def PO2MeanOutputFilePath(self):
        return os.path.join(self.param_study['path'], self.PO2MeanOutputFileName)

    def MTCOutputFilePath(self):
        return os.path.join(self.param_study['path'], self.MTCOutputFileName)

    def probePositionFromIndex(self, probeIdx):
        """ Compute the probe position from the index of the probe."""
        return probeIdx*self.probeSpacing

    def write_results(self):
        with open(self.PO2MeanOutputFilePath(), 'w') as f1,\
             open(self.MTCOutputFilePath(), 'w') as f2:
            po2_writer = csv.writer(f1, delimiter='\t')
            header = list(self.param_study.paramPrefixes())
            header.extend(self.probeNames)
            po2_writer.writerow(header)

            mtc_writer = csv.writer(f2, delimiter='\t')
            header = ['LD', 'U', 'MTC', 'Nu']
            mtc_writer.writerow(header)

            for pp, (LD, U) in itertools.izip(self.get_post_processors(), self.param_study.paramValues()):
                sampled_set = SampledSet(pp.case_path, self.sampledSetDir)
                PO2Values = sampled_set.last_time_values(self.sampledSetName,
                                                        self.sampledField)[:,1]
                stats = ['%g' % LD, '%g' % U]
                stats.extend(PO2Values)
                po2_writer.writerow(stats)

                MTC, Nu = extractMTCEulerian(pp)
                mtc_writer.writerow([LD, U, MTC, Nu])

    def writePO2Mean(self):
        with open(self.PO2MeanOutputFilePath(), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            header = list(self.param_study.paramPrefixes())
            header.extend(self.probeNames)
            writer.writerow(header)
            for pp, (LD, U) in zip(self.get_post_processors(), self.param_study.paramValues()):
                sampled_set = SampledSet(pp.case_path, self.sampledSetDir)
                PO2Values = sampled_set.last_time_values(self.sampledSetName,
                                                        self.sampledField)[:,1]
                stats = ['%g' % LD, '%g' % U]
                stats.extend(PO2Values)
                writer.writerow(stats)

    def writeMTCNu(self):
        with open(self.MTCOutputFilePath(), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            header = ['LD', 'U', 'MTC', 'Nu']
            writer.writerow(header)
            # loop through cases
            for pp, (LD, U) in zip(self.get_post_processors(), self.param_study.paramValues()):
                print 'MTC/Nu for LD = %g and U = %g' % (LD, U)
                MTC, Nu = extractMTCEulerian(pp)
                writer.writerow([LD, U, MTC, Nu])

    def fitIntravascularResistance(self):
        """Fit the intravascular resistance coefficient to simulations.
        """
        data = np.loadtxt(self.PO2MeanOutputFilePath(), delimiter="\t", skiprows=1)
        simulatedPO2 = data[:, self.probeIdx + 2]
        res = minimize_scalar(self._residualFit,
                bounds=[0.1e6, 20e6],
                args=(simulatedPO2),
                method='bounded', options={'disp': True})
        return res.x

    def interpolateTissuePO2(self, LD, U):
        """Interpolate simulated PO2 values to the given values of LD and U.

        Args:
            LD: array-like with values of LD
            U:  array-like with values of U

        Returns:
            numpy array with interpolated values of tissue PO2 to xFit and rFit
        """

        simLD = np.asarray(self.param_study['LD'])
        simU  = np.asarray(self.param_study['U'])
        data = np.loadtxt(self.PO2MeanOutputFilePath(), delimiter="\t", skiprows=1)
        simulatedPO2 = data[:, self.probeIdx + 2]
        simulatedPO2 = np.reshape(np.asarray(simulatedPO2), (len(simLD), len(simU)))
        f = RectBivariateSpline(simLD, simU, simulatedPO2)
        return f(LD, U, grid=False)

    def _residualFit(self, K_IV, simulatedPO2):
        """Compute the residual of the IVR fit using sum of squares."""
        for kroghSol in self.kroghSols:
            kroghSol.intravascularResistanceLDHalf = K_IV
        analyticalTissuePO2 = np.array(map(self._analyticalTissuePO2, self.kroghSols))
        simulatedPO2 = np.asarray(simulatedPO2)
        weights = np.ones(self.param_study.nCases())
        weights[simulatedPO2 < 2] = 0.0
        return np.sum(weights*(simulatedPO2 - analyticalTissuePO2)**2)

    def _analyticalTissuePO2(self, kroghSol):
        return kroghSol.PO2Tissue(self.xFit, self.rFit)

    def _simulatedTissuePO2(self, postProcessor):
        sampled_set = SampledSet(postProcessor.case_path, self.sampledSetDir)
        array = sampled_set.last_time_values(self.sampledSetName, self.sampledField)
        return array[self.probeIdx, 1]
