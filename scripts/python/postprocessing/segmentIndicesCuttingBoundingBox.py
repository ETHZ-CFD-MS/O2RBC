#!/usr/bin/env python
"""
Display the coordinate intervals of the edge segments that cut the domain bounding box.
"""

import argparse

from HbO2.postprocess.case import GraphCasePostProcessor
from utilities.arguments import add_case_argument


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_case_argument(parser)
    args = parser.parse_args()
    case_name = args.caseName
    pp = GraphCasePostProcessor(case_name)
    print 'Segments cutting bounding box:'
    print '\n'.join(['{:d}: '.format(ei) + str(pp.coordinate_intervals_in_box(ei))
                     for ei in pp.rbcDataPostProcessor.edgeIdsWithPaths()])
