#!/usr/bin/env python
#
# Post-process probes: plot EATs and min/max for each probe.

import os
import argparse

# find probe names for PO2
probeNames = []
for f in os.listdir('.'):
    if os.path.isdir(f) and f.startswith('probe') \
                        and f.endswith('PO2'):
        probeNames.append(f)

# probeNames=['UpstreamPO2', 'MidstreamPO2', 'DownstreamPO2']
from_time = 2.0
max_time = 4.0

for probeName in probeNames:
    os.system('plotProbeEATs.py -f PO2 -p %s --from_time %d' %
            (probeName, from_time))
    os.system('plotProbeMinMax.py -f PO2 -p %s --from_time %d' % 
            (probeName, from_time))
    os.system('plot_probe.py -f PO2 -p %s --positions 0 1 2 3 4' \
              '  --max-time %g' % (probeName, max_time))

