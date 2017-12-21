#!/usr/bin/env python
"""
Plot topological and functional radii for different simulations.
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.stats

from HbO2.plot import styles as cosh_styles
from HbO2.plot.tissuevolumes import TissueVolumeCombinedPlotter
from HbO2.postprocess.results import MultiCaseResults
from HbO2.postprocess.tissuevolumes import TissueVolumeResults

from plot.figureoptions import FigureOptions

style_scheme = cosh_styles.StyleSchemeCOSH()
plot_styles = [style_scheme['MVN1'], style_scheme['MVN2']]

parser = argparse.ArgumentParser()
fig_options = FigureOptions(parser)
args = parser.parse_args()

cbf_path = os.environ['OF_CBF']
func_net_path = 'HbO2/eulerGraph/COSH/functionalNetwork'
case_paths_uniform_inlet = [os.path.join(cbf_path, func_net_path,
                                         'funcNet20151213/d_mean_5_d_std_1.5/ROI2b-dx_1um-24proc-new'),
                            os.path.join(cbf_path, func_net_path,
                                         'funcNet20150715/d_mean_5_d_std_1.5/ROI2-dx_1um-24proc-new')]
case_paths_random_inlet = [os.path.join(cbf_path, func_net_path,
                                        'funcNet20151213/d_mean_5_d_std_1.5/ROI2b-dx_1um-24proc-randomInletHb'),
                           os.path.join(cbf_path, func_net_path,
                                        'funcNet20150715/d_mean_5_d_std_1.5/ROI2-dx_1um-24proc-randomInletHb')]

fig_options.parseOptions(args)

results_uniform = MultiCaseResults([TissueVolumeResults(case_path)
                                    for case_path in case_paths_uniform_inlet],
                                    np.hstack)
results_random = MultiCaseResults([TissueVolumeResults(case_path)
                                   for case_path in case_paths_random_inlet],
                                   np.hstack)
plotter = TissueVolumeCombinedPlotter([results_random, results_uniform])

plt.clf()
plotter.plot_histogram_tissue_radii_multipanel(
    functional_colors=[style_scheme['int_func_radii_random_inlet_face_color'],
                       style_scheme['int_func_radii_face_color']],
    functional_labels=['functional (random)', 'functional (constant)']
)
fig_options.saveFig('plotCombinedHistogramTissueRadii')


topol_radii = 1e6*np.array(results_uniform.topological_tissue_radii())
func_uniform_radii = 1e6*np.array(results_uniform.functional_tissue_radii())
func_random_radii = 1e6*np.array(results_random.functional_tissue_radii())
print
print '{:<30s} {:7.5g} +/- {:.5g} um'.format('Topological radii:',
                                               np.mean(topol_radii),
                                               np.std(topol_radii))
print '{:<30s} {:7.5g} +/- {:.5g} um'.format('Functional radii (constant):',
                                               np.mean(func_uniform_radii),
                                               np.std(func_uniform_radii))
print '{:<30s} {:7.5g} +/- {:.5g} um'.format('Functional radii (random):',
                                               np.mean(func_random_radii),
                                               np.std(func_random_radii))

# Statistical tests for the variance of the tissue radii
# F-test
print
print "Statistical test for the variance of the tissue radii"
print "F-test"
n = len(topol_radii)
F = np.var(func_uniform_radii)/np.var(topol_radii)
pvalue = scipy.stats.f.sf(F, n - 1, n - 1)
print "p-value for uniform func. radii vs topol. radii: ", pvalue
F = np.var(func_random_radii)/np.var(topol_radii)
pvalue = scipy.stats.f.sf(F, n - 1, n - 1)
print "p-value for random func. radii vs topol. radii: ", pvalue
F = np.var(func_random_radii)/np.var(func_uniform_radii)
pvalue = scipy.stats.f.sf(F, n - 1, n - 1)
print "p-value for random func. radii vs uniform func. radii: ", pvalue
# Bartlett test
print "Bartlett test"
_, pvalue = scipy.stats.bartlett(func_uniform_radii, topol_radii)
print "p-value for uniform func. radii vs topol. radii: ", pvalue
_, pvalue = scipy.stats.bartlett(func_random_radii, topol_radii)
print "p-value for random func. radii vs topol. radii: ", pvalue
_, pvalue = scipy.stats.bartlett(func_random_radii, func_uniform_radii)
print "p-value for random func. radii vs uniform func. radii: ", pvalue
# Levene test
print "Levene test"
_, pvalue = scipy.stats.levene(func_uniform_radii, topol_radii, center='mean')
print "p-value for uniform func. radii vs topol. radii: ", pvalue
_, pvalue = scipy.stats.levene(func_random_radii, topol_radii, center='mean')
print "p-value for random func. radii vs topol. radii: ", pvalue
_, pvalue = scipy.stats.levene(func_random_radii, func_uniform_radii, center='mean')
print "p-value for random func. radii vs uniform func. radii: ", pvalue
