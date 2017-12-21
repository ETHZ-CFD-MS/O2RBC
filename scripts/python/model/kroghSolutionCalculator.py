#!/usr/bin/env python
"""
Compute Krogh solution at a given location in an axisymmetric region around a capillary.
"""

import argparse

from HbO2.model.kroghSolution import KroghSolution2DCone
from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric


parser = argparse.ArgumentParser()

parser.add_argument('-x', type=float, help='Axial position')
parser.add_argument('-r', type=float, help='Radial position')

args = parser.parse_args()
x = args.x
r = args.r

sim_params = IOHbO2ParametersAxisymmetric('.')
krogh_sol = KroghSolution2DCone(sim_params)

print "Saturation at x:        {:9.8g}".format(krogh_sol.saturationAtX(x))
print "PO2 at (x, r):          {:9.8g}".format(krogh_sol.PO2Tissue(x, r))
