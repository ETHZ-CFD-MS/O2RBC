#!/usr/bin/env python
"""
Plot data from multiple parameter studies for plasma oxygenation.
"""

import argparse
import matplotlib.pyplot as plt
import os

from HbO2.plot.plasmaoxygenation import PlasmaOxygenationParamStudyPlotter
from HbO2.postprocess.factory.parameterstudy import make_param_study_post_processor
from HbO2.setup.parameterStudy import ParameterStudyFactory
from plot.figureoptions import FigureOptions
from plot.utils import annotate_axis_corner

path_to_parameter_studies = "/local/aluecker/OpenFOAM/aluecker-2.3.0/run/cbf/HbO2/eulerAxisymmetric/PO2PlasmaVsHb"
studies = ["paramStudyShort_rc_1.5um_M0_200_2000_LD_0.1_0.9/",
           "paramStudyShort_rc_2.0um_M0_200_2000_LD_0.1_0.9/"]
annotations = ['A', 'B']

parser = argparse.ArgumentParser()
parser.add_argument('--paramFile',
                    help='Path to parameter study file', default='params.json')
parser.add_argument('--manual-labels', action='store_true',
                    help='Whether to set manually the contour labels')
parser.add_argument('--settingsFile', help='Relative path to the file with postprocessing settings',
                    default='postprocess.json')
fig_options = FigureOptions(parser)
args = parser.parse_args()
file_name = args.paramFile
manual_labels = args.manual_labels
settings_file = args.settingsFile
fig_options.parseOptions(args)

fig, axs = plt.subplots(1, 2, sharey=True)
twin_axs = []
for study, ax, text in zip(studies, axs, annotations):
    path_to_study_file = os.path.join(path_to_parameter_studies, study, file_name)
    param_study = ParameterStudyFactory().makeParameterStudy(path_to_study_file)
    postprocessor = make_param_study_post_processor(param_study)
    plotter = PlasmaOxygenationParamStudyPlotter(postprocessor, fig_options)
    plt.sca(ax)
    plotter.contour_delta_po2(manual_labels=manual_labels)
    twin_axs.append(plt.gca())  # hack to get the twin axes
    annotate_axis_corner(ax, text, xy=(0.94, 0.92))

twin_axs[0].set_ylabel('')
twin_axs[0].set_yticklabels([])
axs[1].set_ylabel('')

fig_options.saveFig('plotPlasmaOxygenation_rc_1.5um_2.0um')
