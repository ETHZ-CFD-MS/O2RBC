import argparse
import cPickle as pickle
from itertools import chain

import numpy as np
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

pathToGraphDict      = "domain/constant/graphDict"
pathToEdgeVelocities = "domain/constant/edgeVelocities"
pathToRBCPaths       = "domain/constant/RBCPaths"

class VGMDataTranslator:
    """Class that translates data from the VGM model to OpenFOAM dictionaries that
       can be used as input for the oxygen transport code. These data contain the
       graph definition, RBC velocities on each edge and RBC paths."""

    def __init__(self):
        self.pickleFile = ""
        self.VGMDict = {}
        self.graphDict = {'adjacencyList': [],
                          'vertexIndices': [],
                          'edgeIndices': [],
                          'vertexPositions': []}
        self.edgeVelocitiesDict = {'edgeVelocities': []}
        self.RBCPathsDict = {}

        self.addOptions()

    def addOptions(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('-f', '--file', help='Path to pickle file')
        args = parser.parse_args()
        self.pickleFile = args.file

    def run(self):
        self.readVGMDataFromPickle(self.pickleFile)
        self.extractDictsFromVGMDict()
        self.writeOpenFOAMDictionaries()


    def readVGMDataFromPickle(self,pathToPickle):
        """
        Read data from a pickle file.
        """
        print "Reading pickle file..."
        self.VGMDict = pickle.load(open(pathToPickle, "rb"))

    def extractDictsFromVGMDict(self):
        """
        Extract various dictionaries from VGM data.
        """
        print "Extracting dictionaries from VGM data"
        self.extractGraphDictFromVGMDict()
        self.extractEdgeVelocitiesFromVGMDict()
        self.extractRBCPathsFromVGMDict()

    def extractGraphDictFromVGMDict(self):
        """
        From the dictionary with VGM data, extract the graph data to a dictionary
        whose format can be used to generate OpenFOAM dictionaries.
        """
        # For the adjacency list, transform the tuples to lists
        self.graphDict['adjacencyList'] = [list(edge) for edge in self.VGMDict['adjacencyList']]

        # extract the vertex numbers from the adjacency list
        vertices = list(chain.from_iterable(self.VGMDict['adjacencyList']))
        vertexSet = set(vertices)
        vertices = sorted(list(vertexSet))
        self.graphDict['vertexIndices'] = vertices

        # If a key 'globalEIndex' is present, use the edge indices from the dictionary.
        # Otherwise, the edge indices are defined from 0 to nEdges-1
        if 'edgeIndices' in self.VGMDict:
            self.graphDict['edgeIndices'] = self.VGMDict['edgeIndices'] 
        else:
            self.graphDict['edgeIndices'] = range(len(self.graphDict['adjacencyList']))

        # scale the vertex positions by 1e-6
        self.graphDict['vertexPositions'] = \
                    [[1e-6*x for x in edge] for edge in self.VGMDict['vertexPositions']]

        # copy the edge diameters
        self.graphDict['segmentDiameters'] =  \
                    [[1e-6*x for x in edge] for edge in self.VGMDict['segmentDiameters']]

        # copy the edge points
        if 'edgePoints' in self.VGMDict:
            self.graphDict['edgeType'] = 'polygonal'
            self.graphDict['edgePoints'] = \
                        [[[1e-6*x for x in point] for point in edge] for edge in self.VGMDict['edgePoints']]
        else:
            self.graphDict['edgeType'] = 'straight'

        self.graphDict['tubeOptions'] = {'useEffectiveDiameter': False,
                                         'outerInnerDiameterRatio': 1.25}

    def extractEdgeVelocitiesFromVGMDict(self):
        """
        From the dictionary with VGM data, extract the graph data to a dictionary
        whose format can be used to generate OpenFOAM dictionaries.
        """
        # The velocities produced by VGM have the following structure:
        # it is a list of lists such that
        #   v[timeIndex][edgeIndex]
        # is the velocity on edge "edgeIndex" at time given by time[timeIndex].
        velocityList = self.VGMDict['edgeVelocities']
        signList     = self.VGMDict['sign']
        # velocityList = self.transposeListOfList(velocityList)
        # signList     = self.transposeListOfList(signList)

        self.edgeVelocitiesDict['edgeVelocities'] = \
                [[1e-3*self.VGMDict['time'][i],
                 [1e-3*v*sign for v, sign in zip(velocityList[i], signList[i])]] 
                 for i in range(len(velocityList))]

    def transposeListOfList(self, listOfList):
        array = np.array(listOfList)
        array = np.transpose(array)
        listOfList = array.tolist()
        return listOfList

    def extractRBCPathsFromVGMDict(self):
        """
        From the dictionary with VGM data, extract the graph data to a dictionary
        whose format can be used to generate OpenFOAM dictionaries.
        """
        # Structure of the RBC paths in the VGM dictionary:
        # The RBC paths are entries with keys of the type 1.0, 2.0, 3.0 etc.
        # Each entry contains 3 lists, e.g.
        #   RBCPath[1.0][0] = times
        #   RBCPath[1.0][1] = edge index of RBC at corresponding time
        #   RBCPath[1.0][2] = curvilinear coordinate of RBC at corresponding time
        # The entries are identified by the float type of the keys.
        for key,value in self.VGMDict.iteritems():
            if isinstance(key, int):
                if len(value[0]) > 0:
                    pathName = 'RBCPath%i' % key
                    self.RBCPathsDict[pathName] = {}
                    self.RBCPathsDict[pathName]['index']   = key
                    self.RBCPathsDict[pathName]['times']   = [1e-3*t for t in value[0]]
                    self.RBCPathsDict[pathName]['edges']   = value[1]
                    self.RBCPathsDict[pathName]['sCoords'] = [max(0, 1e-6*s) for s in value[2]]

    def writeOpenFOAMDictionaries(self):
        """
        Write all the OpenFOAM dictionaries with data obtained from VGM.
        """
        print "Writing OpenFOAM dictionaries..."
        self.writeGraphDict()
        self.writeEdgeVelocities()
        self.writeRBCPaths()

    def writeGraphDict(self):
        """
        Write the graph dictionary.
        """
        writeOpenFOAMDictionaryFromDictionary(self.graphDict, pathToGraphDict)

    def writeEdgeVelocities(self):
        """
        Write the edge velocitiy dictionary.
        """
        writeOpenFOAMDictionaryFromDictionary(self.edgeVelocitiesDict, pathToEdgeVelocities)

    def writeRBCPaths(self):
        """
        Write the RBC path dictionary.
        """
        writeOpenFOAMDictionaryFromDictionary(self.RBCPathsDict, pathToRBCPaths)


def writeOpenFOAMDictionaryFromDictionary(dictToWrite, pathToFile):
    """
    Write an OpenFOAM dictionary from a python dictionary
    """
    print "Writing to file", pathToFile
    # open file without reading its body and add dictionary data to it
    dictFile = ParsedParameterFile(pathToFile, noBody=True)
    for key, value in dictToWrite.iteritems():
        dictFile[key] = value
    dictFile.writeFile()
