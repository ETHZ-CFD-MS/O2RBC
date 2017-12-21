#!/usr/bin/env python
"""
Change the end time of a simulation and prepare the case
for restarting at the last written time step.
"""

import argparse
import os
import shutil

from HbO2.setup.case import SimulationParametersFactory
from HbO2.setup.case import copy_const_fields_first_to_last

average_property_filepath = 'fieldAveragingProperties.gz'


parser = argparse.ArgumentParser()
parser.add_argument('-c', '--case', type=str, default='.', help='Case name')
parser.add_argument('-t', '--time-extend', type=float, help='Additional time to simulate')
parser.add_argument('-a', '--averaging-time', type=float, help='Time averaging window length')
args = parser.parse_args()
case_path = args.case
time_extend = args.time_extend
averaging_time = args.averaging_time

sim_params = SimulationParametersFactory().make_sim_params(case_path, use_numerics=True)

domain_path = os.path.join(case_path, 'domain')
last_time_str = copy_const_fields_first_to_last(domain_path)
last_time = float(last_time_str)

# adapt the end time and the averaging start
end_time = last_time + time_extend
if averaging_time > time_extend:
    raise RuntimeError('The averaging window length cannot be smaller than the additional time.')
sim_params['endTime'] = end_time
sim_params['timeStart'] = end_time - averaging_time
sim_params.update_files()
print "endTime:   {:g}".format(sim_params['endTime'])
print "timeStart: {:g}".format(sim_params['timeStart'])

# rename the averaging properties file
try:
    file_path = os.path.join(domain_path, last_time_str, 'uniform', average_property_filepath)
    file_path_bak = file_path + '.bak'
    shutil.move(file_path, file_path_bak)
except:
    print 'Could not move the averaging properties file {:s}'.format(file_path)
