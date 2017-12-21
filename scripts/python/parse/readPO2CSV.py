#!/usr/bin/env python
#
# Read a CSV file with a header line and return the columns in a dictionary
# with corresponding keys
#

import numpy as np
import csv


def readPO2CSV(filePath, delimiter='\t'):
    num_lines = sum(1 for line in open(filePath))
    header_lines = 1
    n_data = num_lines - header_lines
    keys = []

    # read all values into two numpy arrays: time and positions
    with open(filePath, 'r') as csvfile:
        ncol = 0
        i = 0

        reader = csv.reader(csvfile, delimiter=delimiter, quotechar='"')
        for line in reader:
            if i == 0:
                keys = line
                ncol = len(keys)
                data = np.zeros( (n_data, ncol) )
            else:
                data[i-1,:] = [float(x) for x in line]

            i = i+1

    # put results in a dictionary
    results = {}
    for i in range(len(keys)):
        # extract data and compute RBC center
        results[keys[i]] = data[:,i]

    return results

