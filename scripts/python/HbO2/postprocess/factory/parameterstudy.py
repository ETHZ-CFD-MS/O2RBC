"""
Factory for parameter study postprocessors
"""


from HbO2.COSH.postprocess import COSHParamStudyPostProcessor, COSHLDvRBCParamStudyPostProcessor
from HbO2.postprocess.factory.case import load_settings_file
from HbO2.postprocess.parameterstudy import ParameterStudyPostProcessor, LDvRBCParameterStudyPostProcessor
from HbO2.postprocess.plasmaoxygenation import LDjtPlasmaOxygenationParameterStudyPostProcessor
from HbO2.postprocess.po2drop import PO2DropParameterStudyPostProcessor
from HbO2.postprocess.probes import ProbeParameterStudyPostProcessor

default_param_study_pp_class = ParameterStudyPostProcessor
type_to_param_study_pp_class = {'LDvRBC': LDvRBCParameterStudyPostProcessor,
                                'PO2Drop': PO2DropParameterStudyPostProcessor,
                                'plasmaOxygenation': LDjtPlasmaOxygenationParameterStudyPostProcessor,
                                'probes': ProbeParameterStudyPostProcessor,
                                'COSH': COSHParamStudyPostProcessor,
                                'COSHLDvRBC': COSHLDvRBCParamStudyPostProcessor}


def make_param_study_post_processor(param_study, param_file='postprocess.json', **kwargs):
    settings_dict = kwargs.get('settings_dict',
                               load_settings_file(param_study['path'], param_file=param_file))
    post_processor = default_param_study_pp_class(param_study, settings_dict)
    for postproc_dict in settings_dict['postprocessing']:
        if postproc_dict['type'] in type_to_param_study_pp_class:
            post_processor_class = type_to_param_study_pp_class[postproc_dict['type']]
            arg_dict = {key: postproc_dict[key] for key in postproc_dict if key != 'type'}
            post_processor = post_processor_class(post_processor, **arg_dict)
        else:
            print ("Ignoring postprocessing type {:s} while creating postprocessor for "
                   "a parameter study.").format(postproc_dict['type'])
    return post_processor

