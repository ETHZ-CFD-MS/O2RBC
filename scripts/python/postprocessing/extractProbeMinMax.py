#!/usr/bin/env python
#
# Extract the minimum and maximum PO2 at a probe in an Eulerian simulation
#

import argparse

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema

from postprocessing.loadProbesEuler import loadProbes


## Extract the minimum and maximum PO2 at a probe in an Eulerian simulation
# The EATs are returned in a list.
#
# @param caseName Name of the case where probe data are stored
# @param probes dictionary with probe information (produced by loadProbesEuler)
# @param probeIdx index of the probe for minima/maxima are to be computed
# @return Dictionary with min/max information.
def extractProbeMinMax(caseName, probes, probeIdx, from_time=0.0):

    times = probes['times']
    positionFile = 'RBCPositions.txt'
    positionPath = '%s/%s' % (caseName, positionFile)
    probePosition = probes['positions'][probeIdx,:]
    values = probes['values'][:,probeIdx]

    # interpolation object
    f   = interp1d(times, values, kind='linear')

    PO2_min = []
    PO2_max = []
    PO2_mean = []
    LD = []
    EAT = []

    from_idx = np.searchsorted(times, from_time)

    # Extract minima and maxima
    i_max = argrelextrema(values[from_idx:], np.greater)
    i_min = argrelextrema(values[from_idx:], np.less)

    PO2_max = [values[from_idx+i] for i in i_max]
    PO2_min = [values[from_idx+i] for i in i_min]

    return_dict = \
            dict([('PO2_min',  PO2_min),
                  ('PO2_max',  PO2_max)])

    # convert to numpy arrays
    for key in return_dict.keys():
        return_dict[key] = np.asarray(return_dict[key])

    return return_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--field', '-f', help='Field to plot', default='PO2')
    parser.add_argument('--probeName', '-p', help='Name of the probe directory', 
                        default='probeMidstreamPO2')
    parser.add_argument('--from_time', type=float, help='Time from which to plot EATs',
                    default=0.0)

    args            = parser.parse_args()
    fieldName       = args.field
    probeName       = args.probeName
    from_time       = args.from_time

    probes = loadProbes('.', probeName, fieldName)
    minMaxDict = extractProbeMinMax('.', probes, 2, from_time)

    print 'Averaged minimum PO2: %g' % np.mean(minMaxDict['PO2_min'])
    print 'Averaged maximum PO2: %g' % np.mean(minMaxDict['PO2_max'])






