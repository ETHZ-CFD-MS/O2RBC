#!/usr/bin/env python
"""
Plot IVR coefficients.
"""

import argparse
import numpy as np
import os

import matplotlib.pyplot as plt

from HbO2.model.coefficients import intravascularResistanceAnalytical
from HbO2.plot import labels
from HbO2.setup.simulationParameters import IOHbO2SimulationParametersAxisymmetric
from parse.readfile import load_numeric_to_dictionary
from plot.figureoptions import FigureOptions


script_path_env_var = 'OF_SCRIPTS'
radius_rbc_fit_file_path = 'input/IVRFittedOnRadiusRBCNew.txt'
radius_plasma_fit_file_path = 'input/IVRFittedOnRadiusPlasmaNew.txt'
plasma_color = [0, 0, 0.75]

parser = argparse.ArgumentParser()
fig_options = FigureOptions(parser)
args = parser.parse_args()
fig_options.parseOptions(args)

script_path = os.environ[script_path_env_var]
rbc_file_path = os.path.join(script_path, radius_rbc_fit_file_path)
plasma_file_path = os.path.join(script_path, radius_plasma_fit_file_path)
with open(rbc_file_path, 'r') as f:
    data = load_numeric_to_dictionary(f)
    radius_rbc_fit_values = data['radiusRBC']
    ivr_fitted_on_radius_rbc = data['IVR']
with open(plasma_file_path, 'r') as f:
    data = load_numeric_to_dictionary(f)
    radius_plasma_fit_values = data['radiusPlasma']
    ivr_fitted_on_radius_plasma = data['IVR']

sim_params = IOHbO2SimulationParametersAxisymmetric('.')
sim_params['LDMean'] = 0.5
sleeve_thickness = sim_params['radiusPlasma'] - sim_params['radiusRBC']
wall_thickness = sim_params['radiusWall'] - sim_params['radiusPlasma']
analytical_ivr_radius_rbc = np.zeros(ivr_fitted_on_radius_rbc.shape)
analytical_ivr_radius_plasma = np.zeros(ivr_fitted_on_radius_plasma.shape)

for i, rc in enumerate(radius_rbc_fit_values):
    sim_params['radiusRBC'] = rc
    print sim_params['radiusRBC'], sim_params['radiusPlasma']
    analytical_ivr_radius_rbc[i] = intravascularResistanceAnalytical(sim_params)
for i, rp in enumerate(radius_plasma_fit_values):
    sim_params['radiusPlasma'] = rp
    sim_params['radiusRBC'] = rp - sleeve_thickness
    sim_params['radiusWall'] = rp + wall_thickness
    print sim_params['radiusPlasma'], sim_params['radiusRBC'], sim_params['radiusWall']
    analytical_ivr_radius_plasma[i] = intravascularResistanceAnalytical(sim_params)

fig, ax = plt.subplots(1)
ax.plot(1e6*radius_rbc_fit_values, 1e-6*ivr_fitted_on_radius_rbc, 'k.-')
ax.plot(1e6*radius_rbc_fit_values, 1e-6*analytical_ivr_radius_rbc, 'k--')

# plt.gca().set_xlim([min(1e6*radius_rbc_fit_values), max(1e6*radius_rbc_fit_values)])
plt.plot([2.0, 2.0], [0, 10], color='k', linestyle='dotted')
# ax.text(2.03, 4.1, '$r_p = 2.0 \, \muup m$')
ax.set_xlim([min(1e6*radius_rbc_fit_values), 2.5])
ax.set_xticks([1.5, 1.7, 1.9, 2.1, 2.3, 2.5])
labels.setXLabel('r_c', 'um')
labels.setYLabel('K_{IV,0.5}', 'IVR')

axtop = ax.twiny()
plt.plot(1e6*radius_plasma_fit_values, 1e-6*ivr_fitted_on_radius_plasma, '.-', color=plasma_color)
plt.plot(1e6*radius_plasma_fit_values, 1e-6*analytical_ivr_radius_plasma, '--', color=plasma_color)

axtop.set_xlim([min(1e6*radius_plasma_fit_values), max(1e6*radius_plasma_fit_values)])
axtop.set_ylim([2, 7])
labels.setXLabel('r_p', 'um')
axtop.xaxis.label.set_color(plasma_color)
for t1 in axtop.get_xticklabels():
    t1.set_color(plasma_color)

fig_options.saveFig('plotIVR')
