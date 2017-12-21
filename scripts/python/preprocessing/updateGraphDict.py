#! /usr/bin/env python
"""
Convert the key names of a pickle file to the new key names.
"""

import argparse
import cPickle as pickle


key_map = {'diameters': 'segmentDiameters',
           'edgeTuples': 'adjacencyList',
           'vertexCoords': 'vertexPositions'}

parser = argparse.ArgumentParser(description='Update the key names of a pickled graph dictionary.')
parser.add_argument('-g', help='Name of the input pickle', default='graphDict.pkl')
parser.add_argument('-o', help='Name of the output pickle', default='graphDictNew.pkl')
args = parser.parse_args()
old_name = args.g
new_name = args.o

graph_dict = pickle.load(open(old_name, 'rb'))

for old_key, new_key in key_map.iteritems():
    if old_key in graph_dict:
        graph_dict[new_key] = graph_dict.pop(old_key)

pickle.dump(graph_dict, open(new_name, 'wb'))
print "Wrote pickle file with updated keys to {:s}".format(new_name)
