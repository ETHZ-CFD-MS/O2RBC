"""
Interface for OpenFOAM cases
"""

from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory

from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric, \
                                            IOHbO2SimulationParametersAxisymmetric, \
                                            IOHbO2ParametersGraph, \
                                            IOHbO2SimulationParametersGraph, \
                                            IOHbO2ParametersStraightCapillaries, \
                                            IOHbO2SimulationParametersStraightCapillaries
from HbO2.setup.utils import isAxisymmetricCase, isGraphCase, isGraphStraightCapillariesCase


class SimulationParametersFactory(object):
    """
    Construct the subclass of HbO2Parameters that corresponds to a given path.
    """

    def make_sim_params(self, case_path, use_numerics=False):
        if isAxisymmetricCase(case_path):
            if use_numerics:
                return IOHbO2SimulationParametersAxisymmetric(case_path)
            else:
                return IOHbO2ParametersAxisymmetric(case_path)
        elif isGraphCase(case_path):
            if isGraphStraightCapillariesCase(case_path):
                if use_numerics:
                    return IOHbO2SimulationParametersStraightCapillaries(case_path)
                else:
                    return IOHbO2ParametersStraightCapillaries(case_path)
            else:
                if use_numerics:
                    return IOHbO2SimulationParametersGraph(case_path)
                else:
                    return IOHbO2ParametersGraph(case_path)


def copy_const_fields_first_to_last(sName):
    """
    Copy constant fields required for the simulations from the first directory

    Args:
        sName (str): OpenFOAM case

    Returns:
        copied (list), list of copied file names
        dest (str), destination directory
    """
    source=SolutionDirectory(sName,archive=None,paraviewLink=False)

    sDir=source[0]
    dDir=source[-1]

    copied=dDir.copy(sDir,
                     include='*',exclude=[],
                     overwrite=False,
                     mustExist=False,
                     purge=False)
    return copied, dDir.baseName()
