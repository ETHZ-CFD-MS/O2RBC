#!/usr/bin/env python
"""
Plot hemoglobin saturation profiles and related quantities for one or multiple OpenFOAM simulation
in a capillary network.
"""

import argparse
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np

from HbO2.plot.hemoglobingraph import HemoglobinOnSegmentsPlotter
from HbO2.postprocess.factory.case import make_post_processor

from plot.utils import annotate_axis_corner
from plot.figureoptions import FigureOptions
from utilities.arguments import add_case_argument

annotations = [('A', 'B'),
               ('C', 'D'),
               ('E', 'F'),
               ('G', 'H')]

parser = argparse.ArgumentParser()
add_case_argument(parser)
parser.add_argument('--edgeIndices', '-e', type=int, nargs='+', default=[],
                        help='Edge indices where to plot hemoglobin saturation')
parser.add_argument('--nProfileMin', '-n', type=int, default=1,
                    help='Minimal number of profiles to plot')
parser.add_argument('--timeDiffMaxFit', type=float, default=np.inf,
                    help='Maximal time difference used for the linear fit')
figOptions = FigureOptions(parser)
args = parser.parse_args()
case_name = args.caseName
eids = args.edgeIndices
n_profile_min = args.nProfileMin
time_diff_max_fit = args.timeDiffMaxFit
figOptions.parseOptions(args)

postprocessor = make_post_processor(case_name)
plotter = HemoglobinOnSegmentsPlotter(postprocessor, n_profile_min=n_profile_min)

fig, axs = plt.subplots(len(eids), 2)
flattened_axs = [ax for sublist in axs for ax in sublist]
for i, (ei, annot) in enumerate(zip(eids, annotations)):
    plt.sca(axs[i,0])
    plotter.plot_hb_profiles(ei)
    plotter.plot_hb_mean_pm_std(ei)
    plotter.restrict_x_limits_to_defined_values(ei)
    annotate_axis_corner(axs[i,0], annot[0])

    plt.sca(axs[i,1])
    plotter.plot_hb_drop_vs_time_to_previous_rbc(ei, threshold=time_diff_max_fit)
    annotate_axis_corner(axs[i,1], annot[1])

e_string = '_'.join(['{:d}'.format(ei) for ei in eids])
figOptions.saveFig('plotHbProfilesWithDropVsRBCSpacing_e_{:s}'.format(e_string))
