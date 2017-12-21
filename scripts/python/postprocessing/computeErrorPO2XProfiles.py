#!/usr/bin/env python
#
# Compare longitudinal profiles from two simulations and compute the
# error in different norms.
# The first case given as argument is the reference case.
#

import argparse

import numpy as np
from scipy.interpolate import interp1d

from HbO2.parse import readSampleFiles


def computeErrors(caseName1, caseName2, time, setName):

    setDirectory = 'sets_x_axis'

    # load data
    data1 = readSampleFiles.loadSampleFile(caseName1, setDirectory, time, setName)
    data2 = readSampleFiles.loadSampleFile(caseName2, setDirectory, time, setName)

    x_coords = np.asarray(data1[:,0])
    y1       = np.asarray(data1[:,1])
    y2       = np.asarray(data2[:,1])
    print len(x_coords)

    # create interpolation objects
    f1 = interp1d(x_coords, y1, kind='linear', bounds_error=False, fill_value=0.0)
    f2 = interp1d(x_coords, y2, kind='linear', bounds_error=False, fill_value=0.0)
    
    # compute errors
    absMaxError = computeAbsMaxError(f1, f2, x_coords)
    relMaxError = computeRelMaxError(f1, f2, x_coords)
    absL1Error  = computeAbsL1Error(f1, f2, x_coords)
    relL1Error  = computeRelL1Error(f1, f2, x_coords)
    absL2Error = computeAbsL2Error(f1, f2, x_coords)
    relL2Error = computeRelL2Error(f1, f2, x_coords)

    # print to terminal
    print "Absolute error in max-norm = %s" % absMaxError
    print "Relative error in max-norm = %s" % relMaxError
    print "Absolute error in L1-norm  = %s" % absL1Error
    print "Relative error in L1-norm  = %s" % relL1Error
    print "Absolute error in L2-norm  = %s" % absL2Error
    print "Relative error in L2-norm  = %s" % relL2Error

def computeAbsMaxError(f1, f2, x):
    return max(abs(f1(x) - f2(x)))

def computeRelMaxError(f1, f2, x):
    return computeAbsMaxError(f1,f2,x)/max(abs(f1(x)))

def computeAbsL1Error(f1, f2, x):
    # use the trapezoidal rule
    dx = x[1] - x[0]
    return dx*sum(0.5*(abs(f2(x)[:-1]-f1(x)[:-1]) + \
                       abs(f2(x)[1:] -f1(x)[1:] ))) 

def computeRelL1Error(f1, f2, x):
    # use the trapezoidal rule
    dx = x[1] - x[0]
    return computeAbsL1Error(f1,f2,x)/(dx*sum(0.5*(abs(f1(x)[:-1])+abs(f1(x)[1:]))))

def computeAbsL2Error(f1, f2, x):
    # use the trapezoidal rule
    dx = x[1] - x[0]
    return np.sqrt(dx*sum((0.5*((f2(x)[:-1]-f1(x)[:-1])**2 + \
                                (f2(x)[1:] -f1(x)[1:] )**2))))

def computeRelL2Error(f1, f2, x):
    # use the trapezoidal rule
    dx = x[1] - x[0]
    return computeAbsL2Error(f1,f2,x)/   \
           np.sqrt(dx*sum((0.5*(f1(x)[:-1]**2+f1(x)[1:]**2))))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--time', type=float, help='Time for the snapshot')
    parser.add_argument('--set', help='Name of the sample',
                        default='centerline')
    parser.add_argument('--cases', help='Name of the cases to compare', 
                        nargs=2)

    args    = parser.parse_args()
    time    = args.time
    setName = args.set
    cases   = args.cases

    computeErrors(cases[0], cases[1], time, setName)
