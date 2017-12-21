#!/usr/bin/env python
#
# Extract the EATs for each inter-RBC space in a multi-RBC simulation.

import os

import numpy as np
from scipy.interpolate import interp1d

from HbO2.parse import readSampleFiles
from postprocessing import loadPO2Simulation


# extract EATs for all inter-RBC spaces in a multi-RBC simulation for
# a given time.
# The EATs are returned in a list.
def extractEATs(caseName, L_RBC, line_density, time):

    # extract sample data
    setDirectory = 'sets_x_axis'
    setName      = 'centerline'

    sample_data  = readSampleFiles.loadSampleFile(caseName, setDirectory, \
                                                  time, setName)
    x_coords = sample_data[:,0]
    PO2      = sample_data[:,1]

    # extract the intervals between RBCs
    x_intervals = multiRBCGeometry.getInterRBCCenterlineIntervals(L_RBC,
            line_density)

    f = interp1d(x_coords, PO2, kind='cubic')

    # for each interval, compute the EAT
    EATs    = []
    PO2_max = []
    PO2_min = []

    for interval in x_intervals:
        # find the index of first x-value past the RBC front
        idx1 = np.argmax(x_coords >= interval[0])

        # find the index of last x-value before the next RBC rear
        idx2 = np.argmax(x_coords > interval[1]) - 1

        # find the maximum PO2 on the interval
        PO2_max_interior = max(PO2[idx1:idx2])
        PO2_left  = f(interval[0])
        PO2_right = f(interval[1])
        PO2_max.append(max(PO2_max_interior, PO2_left, PO2_right))
        PO2_min.append(min(PO2[idx1:idx2]))

        EATs.append(PO2_max[-1] - PO2_min[-1])

    return EATs, PO2_max, PO2_min

# Return EATs at all time steps for a multi-RBC simulation.
# The returned variable is a dictionary with four entries:
#
# 'time', 'EATs', 'PO2_min', 'PO2_max'
#
# The last three are 2D numpy arrays.
# 1st dimension: time
# 2nd dimension: inter-RBC index
def extractAllEATs(caseName):

    params = loadPO2Simulation.loadMultiRBCParams(caseName)
    line_density = params['line_density']
    L_RBC        = params['L_RBC']

    # extract EATs for each time
    times   = []
    EATs    = []
    PO2_max = []
    PO2_min = []
    for f in os.listdir(caseName + '/sets_x_axis'):
        time = float(f)
        EAT, pmax, pmin = extractEATs(caseName, L_RBC, line_density, time)
        times.append(time)
        EATs.append(EAT)
        PO2_max.append(pmax)
        PO2_min.append(pmin)

    # transform to numpy arrays
    times   = np.asarray(times)
    EATs    = np.asarray(EATs)
    PO2_max = np.asarray(PO2_max)
    PO2_min = np.asarray(PO2_min)

    # sort results according to time
    reorder = times.argsort()
    times = times[reorder]
    EATs = EATs[reorder]
    PO2_max = PO2_max[reorder]
    PO2_min = PO2_min[reorder]
    
    # put the results in a dictionary
    results = {'time': times,
               'EATs': EATs,
               'PO2_max': PO2_max,
               'PO2_min': PO2_min}

    return results

if __name__ == '__main__':
    
    allEATs = extractAllEATs('.')

    for i in xrange(len(allEATs['time'])):
        print 'Time = %0.3f, EATs = ' % allEATs['time'][i],
        print ['%0.2f' % x for x in allEATs['EATs'][i]]


    

