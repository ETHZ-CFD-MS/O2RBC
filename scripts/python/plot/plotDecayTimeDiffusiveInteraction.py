#!/usr/bin/env python
"""
Plot the decay time for diffusive interaction as a function of linear density and
hemoglobin saturation, based the on the linearized model for diffusive interaction.
"""

import argparse
import matplotlib.pyplot as plt

from matplotlib.ticker import LinearLocator, MultipleLocator, FormatStrFormatter
import numpy as np
import os
import string

from HbO2.model.coefficients import DiffusiveInteractionResistanceAnalytical
from HbO2.postprocess.factory.parameterstudy import make_param_study_post_processor
from HbO2.setup.parameterStudy import ParameterStudyFactory
from HbO2.plot import styles
from HbO2.plot.styles import set_COSH_rc_params
from HbO2.plot.labels import setXLabel, setYLabel
from HbO2.plot.diffusiveinteraction import DiffusiveInteractionParameterStudyPlotter
from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric
from plot.figureoptions import FigureOptions
from plot.utils import annotate_axis_corner, contourLevels, add_panel_layout_parser_options

set_COSH_rc_params()

path_to_parameter_studies = "/local/aluecker/OpenFOAM/aluecker-2.3.0/run/cbf/HbO2/eulerGraph/COSH/straightCapillaryArray/diffusiveInteraction"
path_to_sim_theoretical_k_ci = "/local/aluecker/OpenFOAM/aluecker-2.3.0/run/cbf/HbO2/eulerGraph/COSH/straightCapillaryArray/diffusiveInteraction/LDMeanStudy/case_LDMean_0.3"

studies = ["RBCVelocityStudy",
           "LDMeanStudy",
           "M0Study",
           "domainHeightStudy"]
annotations = string.ascii_uppercase

def plot_decay_time(diff_interaction_coeffs, manual_labels=False):
    """
    Plot the decay time for diffusive interaction as a function of linear density
    and hemoglobin saturation

    Args:
        diff_interaction_coeffs (DiffusiveInteractionResistanceAnalytical)
    """
    ld_vals = np.linspace(0.1, 0.9, 81)
    hb_vals = np.linspace(0.1, 0.9, 81)
    xx, yy = np.meshgrid(ld_vals, hb_vals)
    decay_times = diff_interaction_coeffs.decay_time(xx, yy)
    t_levels = contourLevels(decay_times, 0.03, 0.03)
    manual_locations = [(0.221, 0.813), (0.370, 0.700), (0.506, 0.601),
                        (0.648, 0.508), (0.796, 0.397)]
    c_plot = plt.contour(xx, yy, decay_times, t_levels, colors=None, cmap = plt.cm.jet)
    if manual_labels:
        plt.clabel(c_plot, inline=True, inline_spacing=40,
                   fontsize=8, fmt='%g', manual=manual_locations)
    else:
        plt.clabel(c_plot, inline=1, fontsize=7, fmt='%g')
    # print [(text._x, text._y) for text in c_plot.cl]  # display label positions
    setXLabel('LD', '')
    setYLabel('S', '')


parser = argparse.ArgumentParser()
parser.add_argument('--paramFile',
                    help='Path to parameter study file', default='params.json')
parser.add_argument('--manual-labels', action='store_true',
                help='Whether to set manually the contour labels')
parser.add_argument('--modelNames', nargs='+', help='Models for which to plot results',
                    default=['krogh', 'simple'])
add_panel_layout_parser_options(parser)
fig_options = FigureOptions(parser)
args = parser.parse_args()
file_name = args.paramFile
manual_labels = args.manual_labels
model_names = args.modelNames
fig_options.parseOptions(args)

panel_layout = (3, 1) if not args.panelLayout else args.panelLayout
if panel_layout[0] == 1:
    ax1 = plt.subplot(131)
    ax2 = plt.subplot(132)
    ax3 = plt.subplot(133, sharey=ax2)
else:
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)
plt.sca(ax1)
sim_params = IOHbO2ParametersAxisymmetric(path_to_sim_theoretical_k_ci)
diff_interaction_coeffs = DiffusiveInteractionResistanceAnalytical(sim_params)
plot_decay_time(diff_interaction_coeffs, manual_labels)
ax1.yaxis.set_major_locator(MultipleLocator(0.2))

flattened_axs = [ax2, ax2.twiny(), ax3, ax3.twiny()]
style_scheme = styles.StyleSchemeCOSH()
twin_axes_color = style_scheme['twinAxis']['color']
colors = ['k', twin_axes_color, 'k', twin_axes_color]

for study, ax, color in zip(studies, flattened_axs, colors):
    path_to_study_file = os.path.join(path_to_parameter_studies, study, file_name)
    param_study = ParameterStudyFactory().makeParameterStudy(path_to_study_file)
    postprocessor = make_param_study_post_processor(param_study)
    plotter = DiffusiveInteractionParameterStudyPlotter(postprocessor, fig_options, model_names)
    plt.sca(ax)
    plotter.plotComparisonSaturationDifferenceDecayTime(color=color)
    ax.xaxis.label.set_color(color)
    for t1 in ax.get_xticklabels():
        t1.set_color(color)

# adapt axes limits and tick locations
ax2.set_ylim(100, 200)
ax2.yaxis.set_major_locator(MultipleLocator(20))
ax3.set_ylim(100, 200)
ax3.yaxis.set_major_locator(MultipleLocator(20))
ax2.xaxis.set_major_locator(MultipleLocator(0.2))  # RBC velocity
ax3.xaxis.set_major_locator(MultipleLocator(0.4))  # M_0

# if panel_layout[0] == 1:
#     ax3.set_ylabel('')
#     ax3.tick_params(labelleft='off')

annotate_axis_corner(ax1, 'A', xy=(0.90, 0.88))
annotate_axis_corner(ax2, 'B', xy=(0.90, 0.88))
annotate_axis_corner(ax3, 'C', xy=(0.90, 0.88))

fig_options.saveFig('plotDecayTimesDiffusiveInteractionTheoreticalAndTwinAxes')
