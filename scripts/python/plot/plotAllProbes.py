#!/usr/bin/env python
#
# Generate plots for all available probe folders.
# Uses the scripts generateProbePlots.py that calls
# pyFoamTimelinePlot.py.

import os
import shutil
import numpy as np
import argparse
from generateProbePlots import generateProbePlots

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--min-time', type=float, help='Smallest time that should be used', \
                    default = 0.0)
parser.add_argument('--max-time', type=float, help='Largest time that should be used', \
                    default = 0.2)
args = parser.parse_args()

min_time  = args.min_time
max_time  = args.max_time

# get all the probe* folders
probeFolders = []
for f in os.listdir('.'):
    if f.startswith('probe') and f != 'probePlots':
        probeFolders.append(f)

# for each probe folder, generate the corresponding plot

for f in probeFolders:
    # generate the plots
    generateProbePlots(min_time, max_time, f, [])



    



