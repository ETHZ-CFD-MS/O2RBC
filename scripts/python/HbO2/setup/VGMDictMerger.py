import argparse
import cPickle as pickle


class VGMDictMerger:
    """Class that merges a graph dictionary and a RBC data dictionary produced by VGM."""

    def __init__(self):
        self.graphPickle = ""
        self.RBCPickle   = ""
        self.flowPickle  = ""
        self.outputPickle = "VGMDict.pkl"

        self.addOptions()

    def addOptions(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('graph', help='Path to graph pickle')
        parser.add_argument('-R', '--RBC',  default='RBCdict.pkl',  help='Path to RBC pickle')
        parser.add_argument('-f', '--flow', default='flowdict.pkl', help='Path to flow pickle')
        args = parser.parse_args()
        self.graphPickle = args.graph
        self.RBCPickle   = args.RBC
        self.flowPickle  = args.flow

    def run(self):
        print "Loading pickle files..."
        # Contains the graph topology, plus other not-needed (old) data
        graphDict = pickle.load(open(self.graphPickle, 'rb'))

        # Contains the RBC paths
        RBCDict = pickle.load(open(self.RBCPickle, 'rb'))

        # Contains the flow, velocities and signs
        flowDict = pickle.load(open(self.flowPickle, 'rb'))

        # Create a new dictionary
        VGMDict = {'adjacencyList':[]}

        print "Adding data to new dictionary..."

        # Add the graph data
        VGMDict['adjacencyList']    = graphDict['adjacencyList']
        VGMDict['vertexPositions']  = graphDict['vertexPositions']
        VGMDict['segmentDiameters'] = graphDict['segmentDiameters']
        if 'globalEIndex' in graphDict:
            VGMDict['edgeIndices']  = graphDict['globalEIndex']
        else:
            VGMDict['edgeIndices']  = range(len(graphDict['adjacencyList']))

        if 'points' in graphDict:
            VGMDict['edgePoints']   = graphDict['points']

        # Add the flow data
        VGMDict['edgeVelocities']   = flowDict['v']
        VGMDict['time']             = flowDict['time']
        VGMDict['sign']             = flowDict['sign']

        # Add the RBC data
        for key, value in RBCDict.iteritems():
            if isinstance(key, int):
                VGMDict[key] = value

        # Write it to a new pickle file
        print "Writing merged pickle file to %s..." % self.outputPickle
        pickle.dump( VGMDict, open( self.outputPickle, 'wb' ) )

        print "Done"
