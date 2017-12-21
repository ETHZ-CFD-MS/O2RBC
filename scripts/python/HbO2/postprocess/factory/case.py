"""
Factory for postprocessor objects
"""

import json
import os
import warnings

from HbO2.COSH.postprocess import COSHPostProcessor
from HbO2.postprocess.case import CasePostProcessor, GraphCasePostProcessor
from HbO2.postprocess.hemoglobingraph import HemoglobinOnSegmentsPostProcessor, \
                                             HemoglobinOnWholePathsPostProcessor, \
                                             HemoglobinPostProcessorWithIntegrator
from HbO2.postprocess.plasmaoxygenation import PlasmaOxygenationPostProcessor
from HbO2.postprocess.diffusiveinteraction import DiffusiveInteractionPostProcessor
from HbO2.postprocess.po2drop import PO2DropPostProcessor
from HbO2.postprocess.probes import ProbePostProcessor
from HbO2.setup.utils import isAxisymmetricCase

type_to_postprocessor_class = {
    'COSH': COSHPostProcessor,
    'COSHLDvRBC': COSHPostProcessor,
    'diffusiveInteraction': DiffusiveInteractionPostProcessor,
    'hemoglobinOnSegments': HemoglobinOnSegmentsPostProcessor,
    'hemoglobinOnWholePaths': HemoglobinOnWholePathsPostProcessor,
    'hemoglobinWithIntegrator': HemoglobinPostProcessorWithIntegrator,
    'probes': ProbePostProcessor,
    'PO2Drop': PO2DropPostProcessor,
    'plasmaOxygenation': PlasmaOxygenationPostProcessor
}


def make_post_processor(case_path, param_file='postprocess.json', **kwargs):
    """
    Factory for a postprocessor. If not given, loads settings in a dictionary to create a postprocessor
    with various features, The settings dictionary can also be given as a kwarg.

    The postprocessor construction is based on the decorator design pattern.

    Args:
        case_path (str): path to case
        param_file (str, optional): name of the file with postprocessing settings
        **settings_dict (dict, optional): dictionary with postprocessing settings

    Returns:
        CasePostProcessor instance
    """
    settings_dict = kwargs.get('settings_dict',
                               load_settings_file(case_path, param_file=param_file))
    if isAxisymmetricCase(case_path):
        postprocessor = CasePostProcessor(case_path)
    else:
        postprocessor = GraphCasePostProcessor(case_path)
    for postproc_dict in settings_dict['postprocessing']:
        postprocessor_class = type_to_postprocessor_class[postproc_dict['type']]
        arg_dict = {key: postproc_dict[key] for key in postproc_dict if key != 'type'}
        postprocessor = postprocessor_class(postprocessor, **arg_dict)
    return postprocessor


def load_settings_file(case_path, param_file='postprocess.json'):
    """
    Returns a dictionary with postprocessing and plot settings that are read from a json
    file.

    The returned dictionary must have a key 'postprocessing'.

    Args:
        case_path (str): path to case
        param_file (str): name of settings file file

    Returns:
        dict with key 'postprocessing'
    """
    try:
        json_data = open(os.path.join(case_path, param_file))
        data = json.load(json_data)
        json_data.close()
        if 'postprocessing' not in data:
            RuntimeError('No key "postprocessing" in file {:s}'.format(param_file))
        return data
    except IOError:
        warnings.warn('''Returning empty postprocessing dictionary. '''
                      '''Check if the absence of the file {:s} is desired'''.format(param_file))
        return {'postprocessing': []}
