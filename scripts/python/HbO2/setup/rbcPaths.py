import cPickle as pickle
import json
import numpy as np
import os
import random

from HbO2.setup import VGMDataTranslator


class RBCPathGenerator:
    """Generates RBC paths through capillary networks without branchings,
       given a random distribution of RBC linear densities"""

    def __init__(self, case_path, parameterFile='RBCPathParams.json'):
        """
        Constructor that reads RBC path data from a json file.

        In the json file, the entry 'vRBC' can be a float list with number of entries equal
        to the number of edges, or a single float if the RBC velocity is the same in each edge.
        Similarly, the entry 'LDDistribution' can be either a list of dictionaries, or a single
        dictionary if the linear density is the same in each edge.

        Args:
            case_path (str): path to case
            parameterFile (str, optional): name of parameter file
        """
        self.parameterFile = parameterFile
        jsonData=open(os.path.join(case_path, parameterFile))
        data = json.load(jsonData)
        jsonData.close()

        self.RBCLength = data['RBCLength']
        self.LDDistribution = data['LDDistribution']
        self.vRBC = data['vRBC']
        self.graphDictPickle = data['graphDict']
        self.firstTime = data.get('firstTime', 0)
        self.endTime = data['endTime']

        self.graph = pickle.load(open(os.path.join(case_path, self.graphDictPickle), 'rb'))
        n_edges = len(self.graph['adjacencyList'])
        try:
            if len(self.vRBC) != n_edges:
                raise ValueError('Incorrect number of elements for vRBC.')
        except TypeError:  # raised if vRBC does not have __len__ (e.g. if it is a number)
            self.vRBC = [self.vRBC]*n_edges
        if type(data['LDDistribution']) == list:
            if len(data['LDDistribution']) != n_edges:
                raise ValueError('Incorrect number of elements for LDDistribtion.')
        else:
            self.LDDistribution = [self.LDDistribution]*n_edges
        try:
            self.countercurrentEdges = data['countercurrentEdges']
            self.cocurrent = False
        except KeyError:
            self.cocurrent = True
            self.countercurrentEdges = []

    @classmethod
    def fromparser(cls, parser):
        parser.add_argument('--paramFile', help='Path to parameter file',
                            default='RBCPathParams.json')
        args = parser.parse_args()
        return cls('.', args.paramFile)

    def writeAll(self):
        self.writeRBCPaths()
        self.writeEdgeVelocities()

    def writeRBCPaths(self):
        RBCPathsDict = self.generateRBCPaths()
        VGMDataTranslator.writeOpenFOAMDictionaryFromDictionary(
            RBCPathsDict, VGMDataTranslator.pathToRBCPaths)

    def writeEdgeVelocities(self):
        edgeVelocities = self.generateEdgeVelocities()
        edgeVelocitiesDict = {'edgeVelocities': edgeVelocities}
        VGMDataTranslator.writeOpenFOAMDictionaryFromDictionary(
            edgeVelocitiesDict, VGMDataTranslator.pathToEdgeVelocities)

    def generateRBCPaths(self):
        paths = []
        for eI, eLength in enumerate(self.graph['length']):
            nextEntryTime = self.setFirstEntryTime(eI)
            while nextEntryTime < self.endTime:
                paths.append(self.generateRBCPath(nextEntryTime, eI, eLength))
                nextEntryTime += self.drawRBCTimeDifference(eI)
        paths = self.sortRBCPaths(paths)
        RBCPathsDict = {}
        for i, path in enumerate(paths):
            path['index'] = i
            pathName = 'RBCPath%04i' % i
            RBCPathsDict[pathName] = path
        return RBCPathsDict

    def generateRBCPath(self, entryTime, eI, eLength):
        if self.cocurrent or eI not in self.countercurrentEdges:
            sCoords = [0, eLength]
        else:
            sCoords = [eLength, 0]
        edges   = [eI, eI]
        times   = [entryTime, entryTime + eLength/self.vRBC[eI]]
        return {'edges':   edges,
                'sCoords': sCoords,
                'times':   times}

    def sortRBCPaths(self, paths):
        return sorted(paths, key=lambda k: k['times'][0])

    def setFirstEntryTime(self, eI):
        """
        Set the entry time of the first RBC so that RBC do not
        enter at the same time in all capillaries.
        """
        nEdges = len(self.graph['adjacencyList'])
        if nEdges == 4:
            # avoid that edge pairs (0, 3) and (1, 2) have "consecutive" RBCs
            entryOrderMap = {0: 0, 1: 1, 2: 3, 3: 2}
        else:
            entryOrderMap = {i: i for i in range(nEdges)}
        return self.firstTime + entryOrderMap[eI]/float(nEdges) * self.drawRBCTimeDifference(eI)

    def drawRBCTimeDifference(self, eI):
        return self.drawRBCSpacing(eI)/self.vRBC[eI]

    def drawRBCSpacing(self, eI):
        """Return the spacing between RBC centers in the edge eI, based on a given distribution."""
        distribution = self.LDDistribution[eI]
        if distribution['type'] == "constant":
            LDMean = distribution['LDMean']
            return self.spacing(LDMean)
        if distribution['type'] == "discrete":
            LDValues = distribution['LDValues']
            weights  = distribution['weights']
            if len(LDValues) != len(weights):
                ValueError("""The number of elements in LDValues and weights
                    is not the same.""")
            LD = np.random.choice(LDValues, size=1, p=weights)[0]
            return self.spacing(LD)
        if distribution['type'] == "uniformSpacing":
            minLD = distribution['minLD']
            maxLD = distribution['maxLD']
            maxSpacing = self.spacing(minLD)
            minSpacing = self.spacing(maxLD)
            return random.uniform(minSpacing, maxSpacing)
        if distribution['type'] == "logNormalSpacing":
            mu    = distribution['mu']
            sigma = distribution['sigma']
            s     = random.lognormvariate(mu, sigma)
            return (s+1)*self.RBCLength
        else:
            ValueError("Unknown linear density distribution ""%s"""
                       % distribution)

    def generateEdgeVelocities(self):
        times = [0, 100]
        if self.cocurrent:
            velocities = self.vRBC
        else:
            velocities = [-v if eI in self.countercurrentEdges
                          else v for eI, v in enumerate(self.vRBC)]

        return [[time, velocities] for time in times]

    def spacing(self, LD):
        """
        Return the spacing between RBCs for the given linear density.
        """
        return self.RBCLength/LD
