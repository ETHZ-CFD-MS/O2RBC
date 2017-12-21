#!/usr/bin/env python
"""
Reading and writing parameters for HbO2 computations. Some classes interact directly
with OpenFOAM files, such as the controlDict.
"""

import ConfigParser
import copy
import fileinput
import json
import os
import re
import warnings

import numpy as np

from HbO2.setup.geometry import AxisymmetricGeometry, GraphCapillaries
from HbO2.setup.utils import isAxisymmetricCase, isGraphCase
from parse import readfile


class FileInteraction(object):
    """
    Encapsulates to interact with a file that holds parameters.

    This class specifies the path to the file, the suffix of read and write methods
    used by instances of ParameterInterface. It also holds a flag which indicates whether
    data corresponding to the current file was modified or not. Support to require
    or filter keys is provided. Loaded keys are also tracked.

    Attributes:
        path (str): relative path to the file
        method_suffix (str): suffix of the read/write methods employes by ParameterInterface
        modified (bool): flag for modifications of underlying data
        required_keys (list): keys that need to be present in the file
        filter_keys (list): list from which to extract keys in the file.
            If None, all the keys should be extracted.
        loaded_keys (list): list of keys that were loaded from the file.
    """

    def __init__(self, path, method_suffix):
        """
        Construct using path and method suffix information.

        Args:
            path (str): relative path to the file
            method_suffix (str): suffix of the read/write methods employed by ParameterInterface
        """

        self.path = path
        self.method_suffix = method_suffix
        self.modified = False
        self.required_keys = []
        self.filter_keys = None
        self.loaded_keys = []

    def write_method(self):
        return 'write_' + self.method_suffix

    def read_method(self):
        return 'read_' + self.method_suffix

    def are_required_keys_loaded(self):
        """
        Returns where all the required keys are in loaded_keys.
        """
        return all([key in self.loaded_keys for key in self.required_keys])


class ParameterInterface(object):
    """
    General interface between data held by the class and entries in corresponding files.

    The data are held in the dictionary attribute param_dict. When modifying elements
    of this dictionary, it is the client's responsibility to call the method
    set_modified so that the corresponding files are written out. Then, the method
    update_files updates those files that contain modified values.

    Reading and writing with the following methods are supported:
        - key/value extraction using regular expression (function readKeysFromFile)
        - json files

    Attributes:
        path (str): case path
        file_rules (dict): dictionary with parameter file path as key and corresponding
                           ParameterInterface instance as value
        param_dict (dict): dictionary that holds the parameters
    """

    keyValuePattern = readfile.key_space_value_semicolon_pattern

    def __init__(self, case_path, file_rules):
        self.path = case_path
        self.file_rules = file_rules
        self.param_dict = {}
        self.read()

    def set_modified(self, key):
        """Set the modified flag to the file associated to the key"""
        for file, interaction in self.file_rules.iteritems():
            if key in interaction.loaded_keys:
                interaction.modified = True

    def set_not_modified(self):
        """Sets key modification flags to zero."""
        for file, interaction in self.file_rules.iteritems():
            interaction.modified = False

    def read(self):
        for file, interaction in self.file_rules.iteritems():
            interaction.loaded_keys = getattr(self, interaction.read_method())\
                (interaction.path, interaction.filter_keys)
            if not interaction.are_required_keys_loaded():
                raise IOError("""A required key was not found.
                              Required keys: """, interaction.required_keys)

    def read_key_value_dict(self, file_name, filter_keys=None):
        """Read the given keys from the given file and adds the found values
        to self.

        The keys are read using the function readKeysFromFile.

        Args:
            file_name: name of file to read
            filter_keys: iterable with list of keys to extract from. If None,
                         all keys are extracted.
        Returns:
            List of read keys

        """
        with open(os.path.join(self.path, file_name)) as f:
            matches = readfile.readKeysFromFile(f, self.keyValuePattern, filter_keys)
        for key, value in matches.iteritems():
            self.param_dict[key] = value
        return matches.keys()

    def read_json(self, file_name, filter_keys=None):
        """Read the given keys from the given JSON file and adds the found values
        to self.

        Currently, only dictionaries are supported (lists are not).

        Args:
            file_name (str): file name
            filter_keys (list): iterable with the selected keys to extract.
                If None, all keys are extracted.
        """
        fp = open(os.path.join(self.path, file_name))
        jsonData = json.load(fp)
        fp.close()
        loaded_keys = []
        for key in jsonData:
            if filter_keys is None or key in filter_keys:
                loaded_keys.append(key)
                self.param_dict[key] = jsonData[key]
            else:
                print "Skipping unknown key %s" % key
        return loaded_keys

    def read_config_parser(self, file_name, filter_keys=None):
        """Read the given keys from the given file and adds the found values
        to self. This method uses the ConfigParser module.

        Args:
            file_name (str): file name
            filter_keys (list): iterable with the selected keys to extract.
                If None, all keys are extracted.
        """
        config = ConfigParser.RawConfigParser()
        config.read(os.path.join(self.path, file_name))
        loaded_keys = []
        for key in config.options('Main'):
            if filter_keys is None or key in filter_keys:
                loaded_keys.append(key)
                try:
                    self.param_dict[key] = config.getfloat('Main', key)
                except ValueError:
                    self.param_dict[key] = config.get('Main', key)
        return loaded_keys

    def update_files(self):
        """
        Replaces all key values in the corresponding files.

        Only proceeds with replacement if a key corresponding to the file has been
        modified.
        """
        for file, interaction in self.file_rules.iteritems():
            if interaction.modified:
                getattr(self, interaction.write_method()) \
                    (interaction.path, interaction.loaded_keys)

    def write_key_value_dict(self, file_name, keys):
        """
        Write the given keys to the given file using a regular expression.
        """
        filePath = os.path.join(self.path, file_name)
        for line in fileinput.input(filePath, inplace=1, backup=''):
            line = line.rstrip('\n')
            try:
                (key, value) = readfile.readKeyValueFromLine(line,
                                                             self.keyValuePattern)
                if key in keys:
                    pattern = '%s.*;' % key
                    s_after = '{:<20}{:.8g};'.format(key, self.param_dict[key])
                    line = re.sub(pattern, s_after, line)
            except TypeError:
                pass
            print line

    def write_json(self, file_name, keys):
        """
        Write the given keys to the given JSON file.
        """
        jsonData = {}
        for key in keys:
            try:
                jsonData[key] = self.param_dict[key]
            except KeyError:
                pass
        with open(os.path.join(self.path, file_name), 'w') as f:
            json.dump(jsonData, f, sort_keys=True,
                      indent=4, separators=(',', ': '))

    def write_config_parser(self, file_name, keys):
        """
        Write the given keys to the given file using ConfigParser.
        """
        filePath = os.path.join(self.path, file_name)
        with open(filePath, 'w') as f:
            config = ConfigParser.RawConfigParser()
            config.add_section('Main')
            for key in keys:
                if key in self:
                    try:
                        config.set('Main', key, '%g' % self.param_dict[key])
                    except TypeError:
                        config.set('Main', key, str(self.param_dict[key]))
            config.write(f)


class SimulationParameters(object):
    """
    Interface for simulation parameters.

    The class supports update of dependent variables. The dependence relationships are
    held by the dictionary dependentVariables, with the syntax
    dependentVariables[parentVarName] = [childVar1, childVar2, ...]
    """

    dependentVariables = {}

    def __init__(self, param_dict):
        """
        Args:
            param_dict (dict): dictionary with simulation parameter values
        """
        self.param_dict = param_dict
        self.check_cyclic_dependence()
        self.compute_all_dependent_variables()

    def __setitem__(self, key, item):
        """Update variables that are dependent on the modified variable"""
        self.param_dict[key] = item
        self.compute_dependent_variables(key)

    def __getitem__(self, key):
        """Dict-like access to parameter variables"""
        return self.param_dict[key]

    def keys(self):
        """Return dictionary keys"""
        return self.param_dict.keys()

    def compute_all_dependent_variables(self):
        """Compute all the dependent variables."""
        map(self.compute_dependent_variables, self.dependentVariables)

    def compute_dependent_variables(self, parentVar):
        """Compute the variables that depend on the variable parentVar."""
        if parentVar in self.dependentVariables:
            for dependentVar in self.dependentVariables[parentVar]:
                try:
                    getattr(self, 'compute_%s' % dependentVar)()
                except AttributeError:
                    warnings.warn('No modifier found for %s' % dependentVar)
                    pass

    def check_cyclic_dependence(self):
        """
        Checks whether the dependent variables contains any circular dependence.

        Raises:
            ValueError
        """
        visited = {key: False for key in self.dependentVariables}
        for parentVar in self.dependentVariables:
            if self._has_cycle_with_child(visited, parentVar):
                raise ValueError('Cyclic dependence in dependentVariables')

    def _has_cycle_with_child(self, visited, parent):
        """
        Returns true if the parent forms a cycle with one of its children.

        Note that this method is recursive

        Args:
            visited (dict): dictionary holding information whether node was visited
            parent (str): key of the parent to examine

        Returns:
            True if a cycle was found, false otherwise

        """
        if parent in self.dependentVariables:
            for child in self.dependentVariables[parent]:
                if child in visited and visited[child] == True:
                    return True
                else:
                    visited[child] = True
                    cyclic = self._has_cycle_with_child(visited, child)
                    visited[child] = False
                    return cyclic or False
        else:
            return False


class IOSimulationParameters(SimulationParameters):
    """
    Interface for simulation parameter with input/output to files.
    """
    file_rules = {}

    def __init__(self, case_path='.'):
        self.path = case_path
        self.parameter_interface = ParameterInterface(case_path, self.file_rules)
        self.parameter_interface.file_rules = self.file_rules
        super(IOSimulationParameters, self).__init__(self.parameter_interface.param_dict)

    def __setitem__(self, key, item):
        """Track files concerned by the modification and update dependent variables."""
        self.parameter_interface.set_modified(key)
        super(IOSimulationParameters, self).__setitem__(key, item)

    def update_files(self):
        self.parameter_interface.update_files()


class IOHbO2Parameters(IOSimulationParameters):
    """
    Interface for parameters required for HbO2 computations with input/output to
    file initialConditions.
    """
    file_rules = copy.deepcopy(IOSimulationParameters.file_rules)
    file_rules['initialConditions'] = FileInteraction('initialConditions', 'key_value_dict')

    def __init__(self, case_path='.'):
        super(IOHbO2Parameters, self).__init__(case_path)


class IOHbO2SimulationParameters(IOHbO2Parameters):
    """
    Interface for parameters required for HbO2 computations, with support
    for simulation parameters required by OpenFOAM.

    Some OpenFOAM simulation parameters are now specified in a JSON file.
    However, for backwards compatibility, the old ConfigParser is still
    supported. If no JSON file is found, the code reverts to the
    ConfigParser.
    """
    file_rules = copy.deepcopy(IOHbO2Parameters.file_rules)
    file_rules['controlDict'] = FileInteraction('domain/system/controlDict', 'key_value_dict')
    file_rules['controlDict'].filter_keys = ['deltaT', 'endTime', 'timeStart']
    file_rules['simParams'] = FileInteraction('simParams.json', 'json')

    def __init__(self, case_path='.'):
        super(IOHbO2SimulationParameters, self).__init__(case_path)

    def compute_deltaT(self):
        """Compute deltaT from the CFL number and \Delta x."""
        try:
            self['deltaT'] = self['CFL']*self['dx']/self['RBCVelocity']
        except ZeroDivisionError:
            print "The RBCVelocity is zero, not changing deltaT."
        return self['deltaT']


class HbO2ParametersAxisymmetric(SimulationParameters):
    """
    Interface for parameters required for axisymmetric HbO2 computations.
    """

    dependentVariables = {
        'LDMean':       ['MTCBoxXLeft', 'MTCBoxXRight', 'PO2PlasmaInlet'],
        'HbInlet':      ['PO2RBCInlet', 'PO2PlasmaInlet'],
        'domainLength': ['MTCBoxXLeft', 'MTCBoxXRight'],
        'RBCVolume':    ['RBCLength', 'RBCMeshVolume'],
        'radiusRBC':    ['RBCLength'],
        'RBCLength':    ['RBCHalfLength', 'RBCHalfLengthMinus',
                         'MTCBoxXLeft', 'MTCBoxXRight'],
        'radiusPlasma': ['radiusPlasmaMin'],
        'radiusWall':   ['radiusWallMin'],
    }

    def __init__(self, param_dict):
        super(HbO2ParametersAxisymmetric, self).__init__(param_dict)

    def geometry(self):
        """Return a dictionary-like object with the simulation geoemtry."""
        geom = AxisymmetricGeometry()
        for key in geom:
            geom[key] = self[key]
        return geom

    def nRBC(self):
        """
        Compute the number of RBCs from LD, the domain length and the RBC length.
        If the number of RBCs is given in the simulation parameters, use it directly.
        """
        if 'nRBC' in self.keys():
            return int(self['nRBC'])
        else:
            return int(np.ceil(self['LDMean']*self['domainLength']
                               / self['RBCLength'] + 1))

    def compute_RBCLength(self):
        self['RBCLength'] = self['RBCVolume']/(np.pi*self['radiusRBC']**2)

    def compute_RBCHalfLength(self):
        self['RBCHalfLength'] = 0.5*self['RBCLength']

    def compute_RBCHalfLengthMinus(self):
        self['RBCHalfLengthMinus'] = -0.5*self['RBCLength']

    def compute_RBCMeshVolume(self):
        alpha = 5.  # mesh wedge angle in degrees
        self['RBCMeshVolume'] = 5./360*self['RBCVolume']

    def compute_radiusRBC(self):
        self['radiusRBC'] = np.sqrt(self['RBCVolume']/(np.pi*self['RBCLength']))

    def compute_MTCBox(self):
        MTCBoxCenter = 0.5*self['domainLength']
        boxHalfWidth = 0.5*self['RBCLength']/self['LDMean']
        return MTCBoxCenter - boxHalfWidth, MTCBoxCenter + boxHalfWidth

    def compute_MTCBoxXLeft(self):
        self['MTCBoxXLeft']  = self.compute_MTCBox()[0]

    def compute_MTCBoxXRight(self):
        self['MTCBoxXRight'] = self.compute_MTCBox()[1]

    def compute_radiusWallMin(self):
        self['radiusWallMin'] = self['radiusWall'] - 0.01e-6

    def compute_radiusPlasmaMin(self):
        self['radiusPlasmaMin'] = self['radiusPlasma'] - 0.01e-6

    def compute_PO2RBCInlet(self):
        from HbO2.model.kroghSolution import KroghSolution2DCone
        kroghSol = KroghSolution2DCone(self)
        if 'HbInlet' in self.keys():
            self['PO2RBCInlet'] = kroghSol.chem.hillP(self['HbInlet'])
        elif 'PO2RBCInlet' in self and 'HbInlet' not in self:
            self['HbInlet'] = kroghSol.chem.hillS(self['PO2RBCInlet'])

    def compute_PO2PlasmaInlet(self):
        from HbO2.model.kroghSolution import KroghSolution2DCone
        kroghSol = KroghSolution2DCone(self)
        self['PO2PlasmaInlet'] = self['PO2RBCInlet'] \
                               - kroghSol.intravascResistancePO2Drop(0.0)


class IOHbO2ParametersAxisymmetric(IOHbO2SimulationParameters, HbO2ParametersAxisymmetric):
    """
    Interface for parameters required for HbO2 computations in an axisymmetric
    domain with file input/output
    """
    file_rules = copy.deepcopy(IOHbO2SimulationParameters.file_rules)
    file_rules['geometricData'] = FileInteraction('geometricData', 'key_value_dict')
    file_rules['geometricData'].required_keys = \
        ['radiusRBC', 'radiusPlasma', 'radiusWall',
         'radiusTissueLeft', 'radiusTissueRight',
         'domainLength', 'RBCVolume']

    def __init__(self, case_path='.'):
        IOHbO2SimulationParameters.__init__(self, case_path)
        HbO2ParametersAxisymmetric.__init__(self, self.parameter_interface.param_dict)


class IOHbO2SimulationParametersAxisymmetric(IOHbO2ParametersAxisymmetric):
    """
    Interface for parameters required for HbO2 computations in an axisymmetric
    domain, with support for OpenFOAM simulation parameters.
    """
    file_rules = copy.deepcopy(IOHbO2ParametersAxisymmetric.file_rules)

    dependentVariables = copy.deepcopy(IOHbO2ParametersAxisymmetric.dependentVariables)
    dependentVariables['CFL'] = ['deltaT']
    dependentVariables['dx'] = ['deltaT']
    dependentVariables['RBCVelocity'] = ['deltaT']

    def __init__(self, case_path='.'):
        super(IOHbO2SimulationParametersAxisymmetric, self).__init__(case_path)


class HbO2ParametersGraph(SimulationParameters):
    """
    Parameters for oxygen transport simulations in graphs.

    Attributes:
        graph (GraphCapillaries): graph topology
    """
    def __init__(self, sample_dict, graph):
        super(HbO2ParametersGraph, self).__init__(sample_dict)
        self.graph = graph


class IOHbO2ParametersGraph(IOHbO2Parameters):
    """
    Interface for parameters required for HbO2 computations with a graph of capillaries.
    """
    file_rules = copy.deepcopy(IOHbO2Parameters.file_rules)
    file_rules['geometricData'] = FileInteraction('geometricData', 'key_value_dict')
    openfoam_dict = os.path.join('domain', 'constant', 'graphDict')

    def __init__(self, case_path='.'):
        super(IOHbO2ParametersGraph, self).__init__(case_path)
        self.graph = GraphCapillaries.from_openfoam_dict(
            os.path.join(self.path, self.openfoam_dict))

    def sCoordOffset(self):
        """
        Return the offset between the s-coordinates and the x-coordinates.

        More precisely, s = x + offset.

        Returns:
            float, offset
        """
        return 0.0


class IOHbO2SimulationParametersGraph(IOHbO2ParametersGraph, IOHbO2SimulationParameters):
    """
    Interface for parameters required for HbO2 computations in a graph
    with support for OpenFOAM simulation parameters.
    """
    file_rules = copy.deepcopy(IOHbO2ParametersGraph.file_rules)
    file_rules.update(copy.deepcopy(IOHbO2SimulationParameters.file_rules))

    dependentVariables = copy.deepcopy(IOHbO2ParametersGraph.dependentVariables)
    dependentVariables['CFL'] = ['deltaT']
    dependentVariables['dx'] = ['deltaT']
    dependentVariables['RBCVelocity'] = ['deltaT']

    def __init__(self, case_path='.'):
        super(IOHbO2SimulationParametersGraph, self).__init__(case_path)


class IOHbO2ParametersStraightCapillaries(IOHbO2ParametersGraph):
    """
    Interface for parameters required for HbO2 computations in an array of
    straight parallel capillaries.
    """

    dependentVariables = {
        'RBCVolume':         ['RBCLength'],
        'radiusRBC':         ['RBCLength'],
        'radiusTissueLeft':  ['PO2Tissue', 'HbInit'],
        'radiusTissueRight': ['PO2Tissue', 'HbInit'],
        'LDMean':            ['PO2Tissue', 'HbInit'],
        'RBCVelocity':       ['PO2Tissue', 'HbInit'],
        'HbInlet':           ['PO2Tissue', 'HbInit'],
        'PO2Tissue':         ['PO2Plasma'],
        'O2ConsumptionRate': ['PO2Tissue', 'HbInit']
    }

    def __init__(self, case_path):
        super(IOHbO2ParametersStraightCapillaries, self).__init__(case_path)
        try:
            jsonData=open(os.path.join(self.path, 'RBCPathParams.json'))
            data = json.load(jsonData)
            jsonData.close()
            self.cocurrent = ('countercurrentEdges' not in data)
        except IOError:
            self.cocurrent = True

    def geometry(self):
        """Return a dictionary-like object with the simulation geoemtry."""
        geom = AxisymmetricGeometry()
        for key in geom:
            geom[key] = self[key]
        return geom

    def cocurrentFlow(self):
        return self.cocurrent

    def sCoordsInDomain(self):
        """
        Return a tuple with the lower and upper bound of the s-coordinates that are
        within the computational domain.

        It is assumed that the edges are parallel to the x-axis and that the computational
        domain ranges from 0 to domainLength in the x-direction.

        Returns:
            tuple with lower and upper bound
        """
        vertices = self.graph['vertexPositions']
        edgeTuple1 = self.graph['adjacencyList'][0]
        edgeVertices = [vertices[eI] for eI in edgeTuple1]
        minCoord = -edgeVertices[0][0]
        maxCoord = minCoord + self.geometry()['domainLength']
        return minCoord, maxCoord

    def sCoordOffset(self):
        """
        Return the offset between the s-coordinates and the x-coordinates.

        More precisely, s = x + offset.

        Returns:
            offset (float)
        """
        return self.sCoordsInDomain()[0]

    def compute_RBCLength(self):
        self['RBCLength'] = self['RBCVolume']/(np.pi*self['radiusRBC']**2)

    def compute_PO2Tissue(self):
        from HbO2.model.kroghSolution import KroghSolution2DCone
        kroghSol = KroghSolution2DCone(self)
        self['PO2Tissue'] = kroghSol.averagePO2Tissue()

    def compute_PO2Plasma(self):
        self['PO2Plasma'] = self['PO2Tissue']

    def compute_HbInit(self):
        from HbO2.model.kroghSolution import KroghSolution2DCone
        kroghSol = KroghSolution2DCone(self)
        self['HbInit'] = kroghSol.averageSaturation()


class IOHbO2SimulationParametersStraightCapillaries(IOHbO2ParametersStraightCapillaries,
                                                    IOHbO2SimulationParameters):
    """
    Interface for parameters required for HbO2 computations in a straight capillary
    array with support for OpenFOAM simulation parameters.
    """
    file_rules = copy.deepcopy(IOHbO2ParametersStraightCapillaries.file_rules)
    file_rules.update(copy.deepcopy(IOHbO2SimulationParameters.file_rules))

    dependentVariables = copy.deepcopy(IOHbO2ParametersStraightCapillaries.dependentVariables)
    # no dependent variables are defined in HbO2SimulationParameters
    dependentVariables['CFL'] = ['deltaT']
    dependentVariables['dx'] = ['deltaT']
    dependentVariables['RBCVelocity'] = ['deltaT']

    def __init__(self, case_path='.'):
        super(IOHbO2SimulationParametersStraightCapillaries, self).__init__(case_path)


if __name__ == "__main__":
    from HbO2.setup.case import SimulationParametersFactory
    simParams = SimulationParametersFactory().make_sim_params('.', use_numerics=True)
    if isAxisymmetricCase('.'):
        print "Computed deltaT: %g" % simParams['deltaT']
        print simParams
        print "RBC length after setting the radiusRBC to radiusPlasma:"
        simParams['radiusRBC'] = simParams['radiusPlasma']
        print simParams['RBCLength'], simParams['RBCHalfLength']
        print "RBC length after increasing RBC volume:"
        simParams['RBCVolume'] = 69e-18
        print simParams['RBCLength'], simParams['RBCHalfLength']
        simParams.update_files()
    elif isGraphCase('.'):
        print "Computed PO2Tissue: {:g}".format(simParams['PO2Tissue'])
        print "Computed HbInit: {:g}".format(simParams['HbInit'])
        print "Computed deltaT: {:g}".format(simParams['deltaT'])
        simParams.update_files()

