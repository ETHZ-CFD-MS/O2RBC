#!/usr/bin/env python
"""
Compute the hemoglobin saturation along a capillary network using the ODE model.
"""

import argparse
import cPickle as pickle

from HbO2.model.graphintegration import HemoglobinGraphIntegrator

from utilities.arguments import add_case_argument


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_case_argument(parser)
    parser.add_argument('--graphPickle', default='graph.pkl')
    parser.add_argument('--averageHb', action='store_true')
    parser.add_argument('--functionalRadii', action='store_true')
    args = parser.parse_args()
    case_name = args.caseName
    graph_pickle_name = args.graphPickle
    average_hb = args.averageHb
    topological_radii = not args.functionalRadii
    graph = pickle.load(open(graph_pickle_name))
    integrator = HemoglobinGraphIntegrator(case_name, graph, average_hb=average_hb,
                                           topological_radii=topological_radii)
    print "{:<20s}".format("Distal hemoglobin saturation:")
    print '\n'.join([', '.join(['({:8.6g}, {:8.6g})'.format(int_hb.hb, int_hb.weight)
                                for int_hb in hb_distr])
                     for hb_distr in integrator.distal_hb()])
    print "{:<20s}".format("Distal RBC flow: "), integrator.distal_rbc_flow()
    print "{:<20s} {:g}".format("Mean distal hb: ", integrator.average_distal_hb())
    print "{:<20s} {:g}".format("Std distal hb: ", integrator.std_distal_hb())
