#!/usr/bin/env python
"""
Postprocess hemoglobin values from simulations in a graph
"""

import argparse

from HbO2.postprocess.factory.case import make_post_processor
from utilities.arguments import add_case_argument


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_case_argument(parser)
    parser.add_argument('--settingsFile', default='postprocess.json')
    args = parser.parse_args()
    case_name = args.caseName
    settings_file = args.settingsFile
    pp = make_post_processor(case_name, settings_file)
    pp.write_results()
    print 'Wrote results to {:s}'.format(pp.output_file_name)
    print 'Inlet edges: ', pp.rbcDataPostProcessor.rbc_path_analyzer.inlet_edges()
    print 'Outlet edges:', pp.rbcDataPostProcessor.rbc_path_analyzer.outlet_edges()
    print 'Post-converging bifurcation edges:', pp.rbcDataPostProcessor.rbc_path_analyzer.\
                                                post_converging_bifurcation_edges()
