#!/usr/bin/env python
"""
Read RBC positions and return them in a list of dictionaries with
entries 'front', 'hollow' and 'rear'.
"""

import numpy as np

def readRBCPositions(filePath):

    num_lines = sum(1 for line in open(filePath))
    header_lines = 1
    n_times = num_lines - header_lines

    # read all values into two numpy arrays: time and positions
    with open(filePath, 'r') as f:
        ncol = 0
        i = 0
        for line in f:
            if i == 0:
                ncol = len(line.split())
                n_RBC = (ncol-1)/3
                times     = np.zeros( (n_times) )
                positions = np.zeros( (n_times, 3*n_RBC) )
            else:
                values = [float(x) for x in line.split()]
                times[i-1] = values[0]
                positions[i-1,:] = values[1:]

            i = i+1

    # create a list of dictionaries
    results = []
    for i in range(n_RBC):
        # extract data and compute RBC center
        rear   = positions[:,3*i  ]
        hollow = positions[:,3*i+1]
        front  = positions[:,3*i+2]
        center = 0.5*(rear+front)
        results.append({'rear':   rear,
                        'hollow': hollow,
                        'front':  front,
                        'center': center})

    return (results, times)


if __name__ == "__main__":

    readRBCPositions('RBCPositions.txt')

