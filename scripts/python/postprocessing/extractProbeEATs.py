#!/usr/bin/env python
#
# Extract the EATs at a probe in a Eulerian simulation
#

import argparse
import csv

import numpy as np
from scipy.interpolate import interp1d

from postprocessing.loadProbesEuler import loadProbes
from postprocessing.probeUtils import getRBCPassingTimes


## Extract EATs at a probe in an Eulerian RBC simulation.
# The EATs are returned in a list.
#
# @param caseName Name of the case where probe data are stored
# @param probes dictionary with probe information (produced by loadProbesEuler)
# @param centerlineIdx index of the probe for which EATs are to be computed
# @param tissueIdx index of the probe for which EATs are to be computed
# @return Dictionary with EAT information.
def extractEATs(caseName, probes, centerlineIdx, tissueIdx=-1):

    times = probes['times']
    positionFile = 'RBCPositions.txt'
    positionPath = '%s/%s' % (caseName, positionFile)
    probePosition = probes['positions'][centerlineIdx,:]
    centerlineValues = probes['values'][:,centerlineIdx]
    tissueValues     = probes['values'][:,tissueIdx]

    # interpolation object
    f   = interp1d(times, centerlineValues, kind='linear')
    f_t = interp1d(times, tissueValues,     kind='linear')

    # get RBC passing times
    front_times  = getRBCPassingTimes(probePosition, positionPath, 'front')
    hollow_times = getRBCPassingTimes(probePosition, positionPath, 'hollow')
    rear_times   = getRBCPassingTimes(probePosition, positionPath, 'rear')

    # preprocess these arrays: remove useless entries that do not
    # contribute to EATs, so that the first time that we have is a
    # hollow time and the last one is a front time.
    # More precisely, remove the first front/rear time if it is smaller
    # than the first hollow time.
    # Also the last hollow/rear time it they are larger than the last
    # front time.
    if front_times[0] < hollow_times[0]:
        front_times.pop(0)
    if rear_times[0] < hollow_times[0]:
        rear_times.pop(0)
    if rear_times[-1] > front_times[-1]:
        rear_times.pop()
    if hollow_times[-1] > front_times[-1]:
        hollow_times.pop()

    if (len(front_times) != len(hollow_times) or
        len(front_times) != len(rear_times)):
        print len(front_times), len(hollow_times), len(rear_times)
        print 'Warning, the number of passing times for the ' \
              'front/hollow/rear of the RBC are unequal!'
        print 'Was the number of RBCs in the simulation sufficient?'

    # get RBC length in seconds (probe passing time for an entire RBC)
    RBCTime = (rear_times[1] - front_times[0])

    PO2_min = []
    PO2_max = []
    PO2_mean = []
    PO2_rear = []
    PO2_hollow = []
    PO2_front = []
    PO2_tissue = []
    LD = []
    EAT = []

    # print "Number of passing RBCs:", len(rear_times)

    # for all rear times, get the minimum PO2 between the rear time
    # and the next front time
    for i in range(len(rear_times)):
        t_hollow = hollow_times[i]
        t_rear   = rear_times[i]
        t_front  = front_times[i]

        # Interpolation with higher order
        # t_min_interp = max(min(times), t_hollow-0.01)
        # t_max_interp = min(max(times), t_front+0.01)
        # idx_interp = np.where(np.logical_and(times >= t_min_interp, times <= t_max_interp))
        # f2 = interp1d(times[idx_interp], probeValues[idx_interp], kind='cubic')

        # compute the corresponding PO2 values
        PO2_rear.append(f(t_rear))
        PO2_hollow.append(f(t_hollow))
        PO2_front.append(f(t_front))

        # find the minimal PO2 value between t_hollow and t_front
        # tt = np.linspace(t_rear, t_front, 100)
        tt = np.linspace(t_hollow, t_front, 100)
        t_min = tt[np.argmin(f(tt))]
        PO2_min.append(f(t_min))

        # find the corresponding tissue PO2
        PO2_tissue.append(f_t(t_min))

        # find the maximal PO2 value between t_hollow and t_front
        PO2_max.append(max(f(tt)))

        # compute the mean PO2 value between t_hollow and t_front
        PO2_mean.append(np.mean(f(tt)))
        
        # compute EAT
        EAT.append(PO2_max[-1] - PO2_min[-1])

        # find line density
        LD.append(RBCTime/(RBCTime + (t_front - t_rear)))

    return_dict = \
            dict([('t_hollow', hollow_times),
                 ('t_rear', rear_times),
                 ('t_front', front_times),
                 ('PO2_min',  PO2_min),
                 ('PO2_max',  PO2_max),
                 ('PO2_mean',  PO2_mean),
                 ('PO2_rear', PO2_rear),
                 ('PO2_hollow', PO2_hollow),
                 ('PO2_front' , PO2_front),
                 ('PO2_tissue', PO2_tissue),
                 ('LD', LD),
                 ('EAT', EAT)])

    # convert to numpy arrays
    for key in return_dict.keys():
        return_dict[key] = np.asarray(return_dict[key])

    return return_dict


def writeEATToCSV(EAT_dict, probeName):

    outFile = 'EAT%s.csv' % probeName
    fieldnames = ['LD', 'EAT', 'PO2_min', 'PO2_max', 'PO2_hollow', 'PO2_rear',
                  'PO2_front', 't_hollow', 't_rear', 't_front']

    with open(outFile, 'wb') as f:
        writer = csv.DictWriter(f, fieldnames, delimiter='\t')
        # write header
        writer.writerow(dict((fn,fn) for fn in fieldnames))
        # write rows with data
        for i in range(len(EAT_dict['LD'])):
            row = dict((fn, EAT_dict[fn][i]) for fn in fieldnames)
            writer.writerow(row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--field', '-f', help='Field to plot', default='PO2')
    parser.add_argument('--probeName', '-p', help='Name of the probe directory', 
                        default='probeMidstreamPO2')

    args            = parser.parse_args()
    fieldName       = args.field
    probeName       = args.probeName

    probes = loadProbes('.', probeName, fieldName)
    EAT_dict = extractEATs('.', probes, 0)

    writeEATToCSV(EAT_dict, probeName)






