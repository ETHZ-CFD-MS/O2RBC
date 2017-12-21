#!/usr/bin/env python
#
# Functions for extract mass transfer coefficients and Nusselt numbers.
#

import argparse
import operator
import os

import numpy as np

from HbO2.postprocess.case import CasePostProcessor
from postprocessing import readSwakFiles
from postprocessing.probeUtils import getRBCPassingTimes, getClosestRBC

PO2WallSetName = "PO2PlasmaOutFaceCenter"
O2FluxDensitySetName = "O2FluxDensityPlasmaOutCenter"
positionFile = 'RBCPositions.txt'


def extractMTCEulerian(postProcessor):
    simParams = postProcessor.simParams
    L_domain    = simParams['domainLength']
    x_center = 0.5*L_domain
    minSampleTime = 2
    print "Using minSampleTime = {:g}".format(minSampleTime)

    MTCDict = extractMTCAllTimes(postProcessor, is_lagrangian=False)
    O2Flux = MTCDict['flux']
    PO2_RBC_mean = MTCDict['PO2RBCMean']
    PO2_mm = MTCDict['PO2_mm']
    PO2Wall = MTCDict['PO2Wall']

    # select results when a RBC close to the domain center
    probePosition = [x_center, 0, 0]
    MTCDictFiltered = filterByRBCPosition(MTCDict, postProcessor.casePath, probePosition)

    (minIndex, minTime) = closestTimeToCenter(MTCDictFiltered, minSampleTime)
    print "Sample time = ", MTCDictFiltered['times'][minIndex]
    print 'Flux = {:g}, PO2 RBC mean = {:g}, PO2_mm = {:g}, PO2 inner wall = {:g}' \
        .format(O2Flux[minIndex], PO2_RBC_mean[minIndex], PO2_mm[minIndex], PO2Wall[minIndex])

    MTCFinal = MTCDictFiltered['MTC'][minIndex]
    NuFinal  = MTCDictFiltered['Nu'] [minIndex]
    return MTCFinal, NuFinal


def extractMTCLagrangian(postProcessor):
    simParams = postProcessor.simParams
    L_domain    = simParams['domainLength']
    minSampleTime = 0.5
    print "Using minSampleTime = {:g}".format(minSampleTime)

    MTCDict = extractMTCAllTimes(postProcessor, is_lagrangian=True)
    keys = ['times', 'MTC', 'Nu', 'flux', 'PO2RBCMean', 'PO2_mm', 'PO2Wall']
    print '\t'.join(keys)
    print '\n'.join(['\t'.join(['{:g}'.format(MTCDict[key][i]) for key in keys])
                     for i in range(len(MTCDict['times']))])


def extractMTCAllTimes(postProcessor, is_lagrangian = False):
    """
    Compute the mass transfer coefficient and Nusselt number for a giveuxin postprocessor.
    Args:
        postProcessor (casePostProcessor)

    Return:
        tuple with MTC and Nusselt number
    """
    domainPath = os.path.join(postProcessor.case_path, 'domain')
    # read relevant parameters
    simParams = postProcessor.simParams
    L_domain    = simParams['domainLength']
    x_center = 0.5*L_domain
    d_plasma = 2.*simParams['radiusPlasma']
    P50 = simParams['P50']
    hill_exponent = simParams['hillExponent']
    D_p     = simParams['kappaO2Plasma']
    alpha_p = simParams['alphaPlasma']

    # read RBC positions
    RBCPositions = postProcessor.rbc_data.RBCPositions

    # read RBC fields
    RBCFieldDicts = postProcessor.rbc_data.RBCFields

    # load PO2 flux
    filePath = readSwakFiles.createSwakFilePath(domainPath, O2FluxDensitySetName)
    O2FluxArray = readSwakFiles.loadSwakFile(filePath)

    # load averaged PO2
    filePath = readSwakFiles.createSwakFilePath(domainPath, PO2WallSetName)
    PO2WallArray = readSwakFiles.loadSwakFile(filePath)

    times = O2FluxArray[:,0]
    O2Flux = O2FluxArray[:,1]
    PO2Wall = PO2WallArray[:,1]

    PO2_RBC_mean = np.zeros( (len(times)) )
    S_RBC_mean = np.zeros( (len(times)) )
    PO2_mm = np.zeros( (len(times)) )
    MTC = np.zeros( (len(times)) )
    Nu  = np.zeros( (len(times)) )

    for i, time in enumerate(times):
        # get closest RBC
        if is_lagrangian:
            i_closest = 0
        else:
            i_closest = getClosestRBC([x_center], RBCPositions, i)
        # get mean PO2 of that RBC
        PO2_RBC_mean[i] = RBCFieldDicts[i_closest]['PO2_mean'][i]
        # get mixed mean PO2 of that RBC
        S_RBC_mean[i] = RBCFieldDicts[i_closest]['Hb_mean'][i]
        PO2_mm[i] = P50 * pow(S_RBC_mean[i]/(1. - S_RBC_mean[i]), 1./hill_exponent)
        # compute MTC
        MTC[i] = O2Flux[i]/((PO2_RBC_mean[i] - PO2Wall[i]))
        # compute Nusselt number
        Nu[i] = O2Flux[i]*d_plasma/(D_p*alpha_p*(PO2_mm[i] - PO2Wall[i]))
        # print times[i], i_closest, O2Flux[i], PO2_RBC_mean[i], PO2_mm[i], PO2Wall[i], MTC[i], Nu[i]

    return {'MTC': MTC,
            'Nu' : Nu,
            'flux': O2Flux,
            'PO2RBCMean': PO2_RBC_mean,
            'PO2_mm': PO2_mm,
            'PO2Wall': PO2Wall,
            'times': times}


def filterByRBCPosition(MTCDict, caseName, probePosition):
    times = MTCDict['times']
    domainPath = os.path.join(caseName, 'domain')
    positionPath = os.path.join(domainPath, positionFile)
    centerTimes  = getRBCPassingTimes(probePosition, positionPath, 'center')

    filteredIndices = []
    timesToCenter = []
    for centerTime in centerTimes:
        # find the closest time in the array times
        ic = (abs(centerTime - times)).argmin()
        filteredIndices.append(ic)
        timesToCenter.append(centerTime - times[ic])

    MTCDictFiltered = {'MTC':  MTCDict['MTC']   [filteredIndices],
                       'Nu':   MTCDict['Nu']    [filteredIndices],
                       'times':MTCDict['times'] [filteredIndices],
                       'RBCTimesToCenter': timesToCenter}
    return MTCDictFiltered


def closestTimeToCenter(MTCDict, minSampleTime=2.0):
    times         = MTCDict['times']
    timesToCenter = MTCDict['RBCTimesToCenter']
    filteredTimesToCenter = [timesToCenter[i] for i, t in enumerate(times) if t >= minSampleTime]
    filteredIdx           = [i                for i, t in enumerate(times) if t >= minSampleTime]
    minIndex, minTime = min(enumerate(np.abs(filteredTimesToCenter)), key=operator.itemgetter(1))
    minIndex = filteredIdx[minIndex]
    # print "Minimum time to center: ", minTime, timesToCenter[minIndex]
    return (minIndex, minTime)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--lagrangian', help='MTC extraction for Lagrangian simulations',
                        action='store_true')
    args = parser.parse_args()
    is_lagrangian = args.lagrangian
    postProcessor = CasePostProcessor('.')
    if is_lagrangian:
        extractMTCLagrangian(postProcessor)
    else:
        MTC, Nu = extractMTCEulerian(postProcessor)
        print 'MTC = ', MTC
        print 'Nu =  ', Nu
