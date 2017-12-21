"""
Factory for parameter study plotters
"""

from HbO2.COSH.plotCOSH import COSHParamStudyPlotter
from HbO2.plot.diffusiveinteraction import DiffusiveInteractionParameterStudyPlotter
from HbO2.plot.parameterstudy import ParameterStudyPlotter
from HbO2.plot.plasmaoxygenation import PlasmaOxygenationParamStudyPlotter
from HbO2.plot.po2drop import PO2DropParamStudyPlotter
from HbO2.plot.probes import ProbeParamStudyPlotter
from HbO2.postprocess.factory.case import load_settings_file

default_param_study_plotter_class = ParameterStudyPlotter
type_to_param_study_plotter_class = {'COSH': COSHParamStudyPlotter,
                                     'diffusiveInteraction': DiffusiveInteractionParameterStudyPlotter,
                                     'PO2Drop': PO2DropParamStudyPlotter,
                                     'plasmaOxygenation': PlasmaOxygenationParamStudyPlotter,
                                     'probes': ProbeParamStudyPlotter}


def make_param_study_plotter(post_processor, fig_options, param_file='postprocess.json', **kwargs):
    settings_dict = kwargs.get('settings_dict',
                               load_settings_file(post_processor.param_study['path'],
                                                  param_file=param_file))
    try:
        plot_dict = settings_dict['plotting']
        try:
            plotter_class = type_to_param_study_plotter_class[plot_dict['type']]
            arg_dict = {key: plot_dict[key] for key in plot_dict if key != 'type'}
        except KeyError:
            plotter_class = default_param_study_plotter_class
            arg_dict = {}
    except KeyError:
        raise RuntimeError('Missing entry "plotting" in settings dictionary.')
    return plotter_class(post_processor, fig_options, **arg_dict)
