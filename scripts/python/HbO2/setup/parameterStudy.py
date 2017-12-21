"""
Module for the definition of parameter studies.

The relevant parameters are defined in a json file.
"""

import json
import os
import warnings
from abc import ABCMeta, abstractmethod

import numpy as np

from HbO2.model.chemistry import BloodChemistry
from HbO2.setup.case import SimulationParametersFactory
from HbO2.setup.utils import isAxisymmetricCase, isGraphCase


class ParameterStudyFactory(object):
    """
    Class that builds an instance of a derived class of ParameterStudy.

    The class name is read from an input file. If no class name is found,
    the default class name stored as class attribute is used.
    """

    defaultClassName = 'NamedParameterStudy'

    def __init__(self):
        pass

    def makeParameterStudy(self, fileName):
        jsonData = open(fileName)
        data = json.load(jsonData)
        jsonData.close()
        try:
            className = data['className']
        except KeyError:
            className = self.defaultClassName
        try:
            paramStudyClass = globals()[className]
        except KeyError:
            raise KeyError('Class {} not found'.format(className))
        return paramStudyClass(fileName)


class ParameterStudy(dict):
    """
    Base class for parameter studies.

    Implements loading from a json file and construction from a parser.
    """

    __metaclass__ = ABCMeta

    def __init__(self, parameterStudyFile='params.json'):
        super(ParameterStudy, self).__init__()
        self['file'] = parameterStudyFile
        jsonData = open(parameterStudyFile)
        self.data = json.load(jsonData)
        jsonData.close()
        self['path'] = os.path.dirname(self['file'])
        self['baseCaseName'] = self.data['baseCase']
        self['baseCasePath'] = os.path.join(self['path'], self['baseCaseName'])
        self.baseCaseSimParams = SimulationParametersFactory().make_sim_params(self['baseCasePath'])

    @classmethod
    def fromparser(cls, parser):
        parser.add_argument('--paramFile', 
                help='Path to parameter study file', default='params.json')
        args = parser.parse_args()
        return cls(args.paramFile)

    def buildCaseName(self, *args):
        """Construct the case name for given parameter values."""
        prefixes = self.paramPrefixes()
        suffix = '_'.join(['%s_%g' % (prefix, value) for prefix, value 
                                                     in zip(prefixes, args)])
        return 'case_' + suffix

    @abstractmethod
    def paramPrefixes(self):
        """Return a tuple of prefixes for the parameters names."""
        pass

    @abstractmethod
    def paramValues(self):
        """Return the parameter values used in the parameter study.

        Returns:
            a list of tuples with the parameter values.
        """
        return []

    @abstractmethod
    def setParamValues(self, simParams, values):
        """Set the parameters values in simParams based on given parameter values.

        Args:
            simParams: instance of simulationParameter 
            values: tuple of parameter values (should be produced by paramValues(self))
        """
        pass

    def setup(self):
        """Setup case folders for the parameter study.

        Clone the base case and modify the parameters.
        """
        self.cloneCases()
        map(self.setParamValues, self.simParamsList(), self.paramValues())

    def cloneCases(self):
        for casePath in self.casePaths():
            print 'Cloning case %s' % casePath
            os.system('cloneCase.py %s %s' % (self['baseCasePath'], casePath))

    def caseNames(self):
        """Return the list of case names for the parameter study."""
        return [self.buildCaseName(*values) for values in self.paramValues()]

    def casePaths(self):
        """Return the list of case paths for the parameter study."""
        return [os.path.join(self['path'], caseName) for caseName in self.caseNames()]

    def nCases(self):
        """Return the number of cases in the parameter study."""
        return len(self.paramValues())

    def simParamsList(self):
        """Return a list of simulationParameter objects for the parameter study."""
        # return map(lambda s: HbO2SimulationParametersAxisymmetric(s), self.casePaths())
        return map(lambda s: SimulationParametersFactory().\
                             make_sim_params(s, use_numerics=True), self.casePaths())

    def readParameterValues(self, varName):
        """Read the parameter from the parameter study file.

        The values are read from the keys '<varName>Min', '<varName>Max',
        'n<VarName>', or directly from a list of values '<varName>List'.
        """
        try:
            self[varName] = np.linspace(
                                self.data['{}Min'.format(varName)], 
                                self.data['{}Max'.format(varName)], 
                                self.data['n{}{}'.format(varName[0].upper(), 
                                                         varName[1:])])
        except KeyError:
            self[varName] = np.asarray(self.data['{}List'.format(varName)])


class NamedParameterStudy(ParameterStudy):
    """Class for a parameter study constructed from the name of the varied parameter.

    For a parameter named 'varName', an instance reads following keys from the
    json file:
        <varName>Min     
        <varName>Max
        n<VarName>    with capitalized first letter of the variable name.
    """

    def __init__(self, parameterStudyFile='params.json'):
        super(NamedParameterStudy, self).__init__(parameterStudyFile)
        self.varName = self.data['parameterName']
        self.readParameterValues(self.varName)

    def paramPrefixes(self):
        return (self.varName,)

    def paramValues(self):
        return [(x,) for x in self[self.varName]]

    def setParamValues(self, simParams, values):
        simParams[self.varName] = values[0]
        simParams.update_files()


class LDvRBCParameterStudy(ParameterStudy):
    nRBCAverageMinimum = 10  # minimal number of RBCs passing through domain during time averaging

    def __init__(self, parameterStudyFile='params.json'):
        super(LDvRBCParameterStudy, self).__init__(parameterStudyFile)
        self.readParameterValues('LD')
        self.readParameterValues('U')
        self['timeAverageStart'] = self.data['averageStart']

    def paramPrefixes(self):
        return ('LD', 'U')

    def paramValues(self):
        return [(LD, U) for LD in self['LD'] for U in self['U']]

    def setParamValues(self, simParams, values):
        LD = values[0]
        U  = values[1]
        simParams['LDMean'] = LD
        simParams['RBCVelocity'] = U
        timeAverageLength = np.max([0.5, 
                            ((self.nRBCAverageMinimum - 1)*simParams['RBCLength']/LD 
                            + simParams['domainLength'])/U])
        simParams['endTime'] = self['timeAverageStart'] + timeAverageLength
        simParams.update_files()

    def meshGrid(self):
        return np.meshgrid(self['U'], self['LD'])


class LDvRBCParameterStudyConstantMiddlePO2(LDvRBCParameterStudy):
    """
    Parameter study for linear density and RBC velocity that sets the inlet PO2 so that
    PO2 in the middle of the centerline is the same in all simulations, according to the
    analytical model.
    """

    def __init__(self, parameterStudyFile='params.json'):
        super(LDvRBCParameterStudyConstantMiddlePO2, self).__init__(parameterStudyFile)
        self['PO2Middle'] = self.data['PO2Middle']
        self['PO2InletMax'] = self.data['PO2InletMax']

    def setParamValues(self, simParams, values):
        super(LDvRBCParameterStudyConstantMiddlePO2, self).setParamValues(simParams, values)
        from HbO2.model.kroghSolution import KroghSolution2DCone
        kroghSol = KroghSolution2DCone(simParams)
        xMiddle = 0.5*simParams.geometry()['domainLength']
        hbMiddle = kroghSol.chem.hillS(self['PO2Middle'])
        hbInlet = kroghSol.saturationAtXWithO2ConvectionAndInitialCondition(xMiddle, hbMiddle, 0)
        PO2RBCInlet = kroghSol.chem.hillP(hbInlet)
        if PO2RBCInlet > self['PO2InletMax']:
            warnings.warn('RBC PO2 at the inlet is clamped to {:g}'.format(self['PO2InletMax']))
        simParams['PO2RBCInlet'] = min(PO2RBCInlet, self['PO2InletMax'])
        simParams.update_files()


class LDM0ParameterStudy(ParameterStudy):
    nRBCAverageMinimum = 10  # minimal number of RBCs passing through domain during time averaging

    def __init__(self, parameterStudyFile='params.json'):
        super(LDM0ParameterStudy, self).__init__(parameterStudyFile)
        self.readParameterValues('LD')
        self.readParameterValues('O2ConsumptionRate')
        self['timeAverageStart'] = self.data['averageStart']

    def paramPrefixes(self):
        return ('LD', 'O2ConsumptionRate')

    def paramValues(self):
        return [(LD, M0) for LD in self['LD'] for M0 in self['O2ConsumptionRate']]

    def setParamValues(self, simParams, values):
        LD = values[0]
        M0 = values[1]
        simParams['LDMean'] = LD
        simParams['O2ConsumptionRate'] = M0
        U = simParams['RBCVelocity']
        timeAverageLength = np.max([0.5,
                                    ((self.nRBCAverageMinimum - 1)*simParams['RBCLength']/LD
                                     + simParams['domainLength'])/U])
        simParams['endTime'] = self['timeAverageStart'] + timeAverageLength
        simParams.update_files()

    def meshGrid(self):
        return np.meshgrid(self['O2ConsumptionRate'], self['LD'])


class LDM0ParameterStudyConstantMiddlePO2(LDM0ParameterStudy):
    """
    Parameter study for linear density and oxygen consumption rate that sets the inlet PO2
    so that PO2 in the middle of the centerline is the same in all simulations, according to the
    analytical model.
    """

    def __init__(self, parameterStudyFile='params.json'):
        super(LDM0ParameterStudyConstantMiddlePO2, self).__init__(parameterStudyFile)
        self['PO2Middle'] = self.data['PO2Middle']
        self['PO2InletMax'] = self.data['PO2InletMax']

    def setParamValues(self, simParams, values):
        super(LDM0ParameterStudyConstantMiddlePO2, self).setParamValues(simParams, values)
        from HbO2.model.kroghSolution import KroghSolution2DCone
        kroghSol = KroghSolution2DCone(simParams)
        xMiddle = 0.5*simParams.geometry()['domainLength']
        hbMiddle = kroghSol.chem.hillS(self['PO2Middle'])
        hbInlet = kroghSol.saturationAtXWithO2ConvectionAndInitialCondition(xMiddle, hbMiddle, 0)
        PO2RBCInlet = kroghSol.chem.hillP(hbInlet)
        if PO2RBCInlet > self['PO2InletMax']:
            warnings.warn('RBC PO2 at the inlet is clamped to {:g}'.format(self['PO2InletMax']))
        simParams['PO2RBCInlet'] = min(PO2RBCInlet, self['PO2InletMax'])
        simParams.update_files()


class CapillaryRadiusParameterStudy(NamedParameterStudy):

    varName = 'radiusPlasma'

    def __init__(self, parameterStudyFile='params.json'):
        ParameterStudy.__init__(self, parameterStudyFile)
        self.readParameterValues(self.varName)
        self.plasmaSleeveThickness = self.data['plasmaSleeveThickness']
        self.wallThickness         = self.data['wallThickness']

    def setParamValues(self, simParams, values):
        radiusPlasma = values[0]
        simParams['radiusRBC']  = radiusPlasma - self.plasmaSleeveThickness
        simParams['radiusWall'] = radiusPlasma + self.wallThickness
        super(CapillaryRadiusParameterStudy, self).setParamValues(simParams, values)


class SigmaSParameterStudy(NamedParameterStudy):

    def __init__(self, parameterStudyFile='params.json'):
        super(SigmaSParameterStudy, self).__init__(parameterStudyFile)
        self['SMeanInlet'] = self.data['SMeanInlet']

    def paramPrefixes(self):
        return ('sigmaS',)

    def paramValues(self):
        return [(sigmaS,) for sigmaS in self['sigmaS']]

    def setParamValues(self, simParams, values):
        sigmaS = values[0]
        chem = BloodChemistry()
        simParams['PO2RBCInlet']  = chem.hillP(self['SMeanInlet'])
        simParams['PO2InletHigh'] = chem.hillP(self['SMeanInlet'] + sigmaS)
        simParams['PO2InletLow']  = chem.hillP(self['SMeanInlet'] - sigmaS)
        simParams.update_files()


class DeltaSParameterStudyGraph(NamedParameterStudy):

    def __init__(self, parameterStudyFile='params.json'):
        super(DeltaSParameterStudyGraph, self).__init__(parameterStudyFile)
        self['HbMeanInlet'] = self.data['HbMeanInlet']

    def paramPrefixes(self):
        return ('deltaS',)

    def paramValues(self):
        return [(deltaS,) for deltaS in self['deltaS']]

    def setParamValues(self, simParams, values):
        deltaS = values[0]
        simParams['HbInlet1'] = self['HbMeanInlet'] + deltaS/2
        simParams['HbInlet2'] = self['HbMeanInlet'] - deltaS/2
        simParams.update_files()


class CapillarySpacingParameterStudy(NamedParameterStudy):

    def __init__(self, parameterStudyFile='params.json'):
        super(CapillarySpacingParameterStudy, self).__init__(parameterStudyFile)

    def paramPrefixes(self):
        return ('H',)

    def paramValues(self):
        return [(h,) for h in self['H']]

    def setParamValues(self, simParams, values):
        H = values[0]
        simParams['radiusTissueLeft'] = H/np.sqrt(np.pi)
        simParams['radiusTissueRight'] = simParams['radiusTissueLeft']
        simParams.update_files()
