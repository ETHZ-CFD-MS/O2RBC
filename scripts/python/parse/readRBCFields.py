#!/usr/bin/env python
"""
Functions for reading RBC fields in folder RBCData.
"""

import os
import numpy as np
import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def readRBCFields(filePath):
    """ Read RBC fields stored in file 'filePath' and return them in a
       dictionary with entries 'time', 'PO2_mean', 'PO2_max', 'PO2_min',
       'Hb_mean', 'Hb_max', 'Hb_min'."""

    num_lines = sum(1 for line in open(filePath))
    header_lines = 1
    n_times = num_lines - header_lines

    times  = np.zeros( (n_times) )
    fields = np.zeros( (n_times, 6) )

    # read all values into two numpy arrays: time and positions
    with open(filePath, 'r') as f:
        ncol = 0
        i = 0
        for line in f:
            if i > 0:
                values = [float(x) for x in line.split()]
                times[i-1]    = values[0]
                fields[i-1,:] = values[1:]

            i = i+1

    return {'times'   : times,
            'Hb_mean' : fields[:,0],
            'Hb_min'  : fields[:,1],
            'Hb_max'  : fields[:,2],
            'PO2_mean': fields[:,3],
            'PO2_min' : fields[:,4],
            'PO2_max' : fields[:,5]}

def readAllRBCFields(RBCDataPath):
    """ Read fields for all RBCs and return a list of dictionaries produced
        by the function readRBCFields.
        RBCDataPata = path to folder domain/RBCData."""

    RBCFileNames = []
    for f in os.listdir(RBCDataPath):
        if f.startswith('RBC') and \
           f.endswith('.txt'):
            RBCFileNames.append(f)
    
    # sort the file names using natural sorting
    RBCFileNames.sort(key=natural_keys)

    RBCFieldDicts = []

    for i in range(len(RBCFileNames)):
        filePath = os.path.join(RBCDataPath, RBCFileNames[i])
        RBCFieldDicts.append(readRBCFields(filePath))

    return RBCFieldDicts
