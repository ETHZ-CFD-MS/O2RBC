#!/usr/bin/env python
"""
Generate a "graph" with a square array of parallel straight capillaries
with a given length.
"""

import argparse
import cPickle as pickle

import numpy as np

from HbO2.setup import VGMDataTranslator

parser = argparse.ArgumentParser()
parser.add_argument('-L', '--length', type=float, help='Domain length')
parser.add_argument('-H', '--height', type=float, help='Domain height and depth')
parser.add_argument('-d', '--diameter', type=float, default=4e-6, help='Vessel diameter')
parser.add_argument('-w', '--wall-ratio', type=float, default=1.25,
                    help='Diameter ratio between capillary wall and lumen')
parser.add_argument('-n', type=int,   default=2, help='Number of capillaries per side')
parser.add_argument('--capillaryAtOrigin', action='store_true',
                    help='Whether one of the capillary should be at y=0, z=0.')
parser.add_argument('--symmetricDomain', action='store_true',
                    help='Whether a symmetric domain should be used')
parser.add_argument('-o', '--output', default='graphDict.pkl', help='Output pickle file')
args = parser.parse_args()
L = args.length
H = args.height
diameter = args.diameter
wall_diameter_ratio = args.wall_ratio
nSide = args.n
capillaryAtOrigin = args.capillaryAtOrigin
symmetricDomain = args.symmetricDomain
graphDictPickle = args.output

if nSide <= 0:
    ValueError("Invalid value %i for n, should be a positive integer" % nSide)

nEdges = nSide**2
nSide = int(np.sqrt(nEdges))

def vertexPositions(nSide, L, H):
    x0 = -0.1*L
    x1 =  1.1*L
    if capillaryAtOrigin:
        if symmetricDomain:
            yLocations  = [2*H*i/nSide for i in range(nSide)]
        else:
            yLocations  = [H*i/nSide for i in range(nSide)]
    else:
        if symmetricDomain:
            yLocations  = [(1+2*i)/(2.*nSide)*2*H for i in range(nSide)]
        else:
            yLocations  = [(1+2*i)/(2.*nSide)*H for i in range(nSide)]
    zLocations  = yLocations
    return [[x, y, z] for z in zLocations for y in yLocations for x in [x0, x1]] 

def adjacencyList(nEdges):
    return [[2*i, 2*i+1] for i in range(nEdges)]

def edgeLengths(graphDict):
    positions = np.asarray(graphDict['vertexPositions'])
    return [np.linalg.norm(positions[e[1],:] - positions[e[0],:]) 
                for e in graphDict['adjacencyList']]

graphDict = {}
graphDict['vertexPositions'] = vertexPositions(nSide, L, H)
graphDict['vertexIndices'] = range(2*nEdges)
graphDict['edgeIndices'] = range(nEdges)
graphDict['adjacencyList'] = adjacencyList(nEdges)
graphDict['segmentDiameters'] = [[diameter]]*nEdges
graphDict['edgeType'] = 'straight'
graphDict['tubeOptions'] = {'useEffectiveDiameter': False,
                            'outerInnerDiameterRatio': wall_diameter_ratio}
graphDict['length'] = edgeLengths(graphDict)

pickle.dump(graphDict, open(graphDictPickle, 'wb'))
print "Wrote pickle file %s" % graphDictPickle

VGMDataTranslator.writeOpenFOAMDictionaryFromDictionary(graphDict, VGMDataTranslator.pathToGraphDict)
print "Wrote file %s" % VGMDataTranslator.pathToGraphDict
