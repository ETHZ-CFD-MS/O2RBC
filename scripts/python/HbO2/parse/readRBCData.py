"""
Functions for reading RBC data produced by the OpenFOAM code
for axisymmetric geometries.
"""

import os
import re
import warnings

import numpy as np

from parse.readfile import load_csv_to_dictionary


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    """
    return [ atoi(c) for c in re.split('(\d+)', text) ]


def readRBCFields(filePath, minTime=0.0, includeMinMax=False):
    """
    Read RBC fields stored in file 'filePath' and return them in a
    dictionary with entries 'time', 'PO2_mean', 'PO2_max', 'PO2_min',
    'Hb_mean', 'Hb_max', 'Hb_min'.
    """
    with open(filePath, 'r') as f:
        data = load_csv_to_dictionary(f)
    timeIndices = np.where(data['//- Time'] >= minTime)[0]
    results = {'times'   : data['//- Time'][timeIndices],
               'Hb_mean' : data['Mean Hb'][timeIndices],
               'PO2_mean': data['Mean PO2'][timeIndices]}
    if includeMinMax:
        results['Hb_min'] = data['Min Hb'][timeIndices],
        results['Hb_max'] = data['Max Hb'][timeIndices],
        results['PO2_min'] = data['Min PO2'][timeIndices],
        results['PO2_max'] = data['Max PO2'][timeIndices],
    return results


def readAllRBCFields(RBCDataPath, minTime=0.0):
    """
    Read fields for all RBCs into a list of dictionaries

    Args:
        RBCDataPath: path to folder domain/RBCData.

    Returns:
        list of dictionaries produced by the function readRBCFields.
    """
    RBCFileNames = [f for f in os.listdir(RBCDataPath)
                    if f.startswith('RBC') and f.endswith('.txt')]
    if not RBCFileNames:
        warnings.warn('No RBC data files found in {}'.format(RBCDataPath))
    RBCFileNames.sort(key=natural_keys)
    RBCFieldDicts = []
    for fileName in RBCFileNames:
        filePath = os.path.join(RBCDataPath, fileName)
        RBCFieldDicts.append(readRBCFields(filePath, minTime=minTime))
    return RBCFieldDicts


def readRBCPositions(filePath, minTime=0.0):
    """
    Read the positions of RBCs from the tab-separated file with the structure of RBCPositions.txt.

    Args:
        filePath (str): path to file with RBC positions

    Returns:
        tuple (results, times), where results is a list of dictionaries
        with keys 'rear', 'hollow', 'front' and 'center', and times is an array
        with the corresponding time values.
    """
    with open(filePath, 'r') as f:
        data = load_csv_to_dictionary(f)
    nRBC = (len(data.keys()) - 1)/3
    timeIndices = np.where(data['//- Time'] >= minTime)[0]
    times = data['//- Time'][timeIndices]

    results = []
    for i in range(nRBC):
        positionDict = {}
        # positionDict['rear'] = data['Rear{:d}'.format(i)]
        # positionDict['hollow'] = data['Hollow{:d}'.format(i)]
        # positionDict['front'] = data['Front{:d}'.format(i)]
        # positionDict['center'] = 0.5*(positionDict['front'] + positionDict['rear'])
        positionDict['center'] = 0.5*(data['Rear{:d}'.format(i)] + data['Front{:d}'.format(i)])
        positionDict['center'] = positionDict['center'][timeIndices]
        results.append(positionDict)
    return (results, times)
