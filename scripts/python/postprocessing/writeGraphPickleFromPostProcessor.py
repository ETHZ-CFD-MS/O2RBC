#!/usr/bin/env python
"""
Compute the hemoglobin saturation along a capillary network using the ODE model.
"""

import argparse
import cPickle as pickle

from HbO2.setup.graph import make_graph_from_post_processor, plot_graph
from HbO2.postprocess.factory.case import make_post_processor

from utilities.arguments import add_case_argument


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_case_argument(parser)
    parser.add_argument('--settingsFile', default='postprocess.json')
    parser.add_argument('--pickleName', default='graph.pkl')
    args = parser.parse_args()
    case_name = args.caseName
    param_file = args.settingsFile
    pickle_name = args.pickleName
    pp = make_post_processor(case_name, param_file)
    graph = make_graph_from_post_processor(pp)
    plot_graph(graph)
    pickle.dump(graph, open(pickle_name, 'wb'))
    print "Wrote graph pickle to {:s}".format(pickle_name)
