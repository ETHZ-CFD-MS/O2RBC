#!/usr/bin/env python
#
# Load all x-profiles in a multi-RBC simulation.

import numpy as np
import os

from HbO2.parse import readSampleFiles
from postprocessing import loadPO2Simulation

sampleDirectory = 'sets_x_axis'


# return sample values for a given sample at a given time
def loadProfile(caseName, time, sampleName):

    # extract sample data

    sample_data  = readSampleFiles.loadSampleFile(caseName, sampleDirectory,
                                                  time, sampleName)

    return sample_data

# return all sample values, for all samples in sets_x_axis and all
# times. The results are returned in a dictionary with keys
#
# 'times', 'x_coord' and '<profName>',
#
# where <profName> corresponds to any x-profile name. These entries
# are 2D numpy array index by [time, x]
def loadAllProfiles(caseName):

    # get all output times
    times = loadPO2Simulation.returnOutputTimes(caseName)

    # get all x-profile names
    sampleNames = []
    for f in os.listdir(caseName + '/' + sampleDirectory + '/0'):
        # remove the extension '_PO2.xy'
        sampleName = f.replace('_PO2.xy', '')
        sampleNames.append(sampleName)

    # put the results into a dictionary
    results = {'times': times}
    for sampleName in sampleNames:
        profiles = []
        for t in times:
            data = loadProfile(caseName, t, sampleNames[0])
            x_coord = data[:,0]
            profiles.append(data[:,1])

        profiles = np.asarray(profiles)
        results[sampleName] = profiles

    results['x_coord'] = x_coord

    return results

