#!/usr/bin/env python
#
# Extract the statistical moments of values at a probe.
#

import argparse
import numpy as np

from postprocessing.loadProbesEuler import loadProbes

## Extract the statistical moments of values at a probe.
#
# @param caseName Name of the case where probe data are stored
# @param probes dictionary with probe information (produced by loadProbesEuler)
# @param probeIdx index of the probe for which statistics are to be computed
# @return Dictionary with min/max information.
def extractProbeMoments(caseName, probes, probeIdx, fromTime=0.0):

    times = probes['times']
    values = probes['values'][:,probeIdx]
    nProbes = probes['nProbes']

    sampleTimeIdx = [i for i,t in enumerate(times) if t >= fromTime]

    mean     = [np.mean(probes['values'][sampleTimeIdx, i]) for i in range(nProbes)]
    std      = [np.std (probes['values'][sampleTimeIdx, i]) for i in range(nProbes)]
    variance = [np.var (probes['values'][sampleTimeIdx, i]) for i in range(nProbes)]

    return_dict = \
            dict([('mean', mean    ),
                  ('std',  std     ),
                  ('var',  variance)])

    # convert to numpy arrays
    for key in return_dict.keys():
        return_dict[key] = np.asarray(return_dict[key])

    return return_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--field', '-f', help='Field to plot', default='PO2')
    parser.add_argument('--probeName', '-p', help='Name of the probe directory')
    parser.add_argument('--fromTime', type=float, help='Time from which extract statistics', 
                    default=0.0)

    args            = parser.parse_args()
    fieldName       = args.field
    probeName       = args.probeName
    fromTime        = args.fromTime

    probes = loadProbes('.', probeName, fieldName)

    probePositions = range(probes['nProbes'])

    statDict = extractProbeMoments('.', probes, probePositions, fromTime)


    for i in range(probes['nProbes']):
        positionStr = '(%5.2e, %5.2e)' % (probes['positions'][0,i],
                                          probes['positions'][1,i])
        print 'Probe at position %s:' % positionStr
        print 'Mean: %g' % statDict['mean'][i]
        print 'Std : %g' % statDict['std'][i]
        print 'Var : %g' % statDict['var'][i]







