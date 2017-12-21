#!/usr/bin/env python
#
# Reads a dictionary and adds some stuff to it.
#
# Usage: writeToDictionary.py
#

import numpy as np

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

dictPath='sandboxDict'

vertexIndices = range(0,4)
edgeIndices = range(0,3)
edgeTuples=[[1,2],[2,3],[2,4]]

vertexPositions = [[1,0,0],[2,1,1],[3,2,2],[0,0,-1.5]]
edgeDiameters = [4] * 3
edgeFlows = [1] * 3

RBCPath1 = {};
RBCPath1['index'] = 1
RBCPath1['times'] = [0,0.5,1]
RBCPath1['edges'] = [0,0,1]
RBCPath1['sCoords'] = [0.25,0.33,0.15]

file = ParsedParameterFile(dictPath)


# assert file['model']=='theModel'
# assert len(file['listValues']) > 1
# print "Assertions successfully passed"

file['subDict'] = {"a":1, "b": 'asd'}

file['vertexIndices'] = vertexIndices
file['edgeIndices'] = edgeIndices
file['adjacencyList'] = edgeTuples

file['vertexPositions'] = vertexPositions
file['edgeDiameters'] = edgeDiameters
file['edgeFlows'] = edgeFlows

file['RBCPath1'] = RBCPath1


print file


print "End"
