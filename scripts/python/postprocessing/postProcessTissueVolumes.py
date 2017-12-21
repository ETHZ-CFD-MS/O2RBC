#!/usr/bin/env python
"""
Postprocess tissue volumes obtained from simulations in a graph
"""

import argparse
import numpy as np
from scipy.stats.stats import linregress

from HbO2.postprocess.results import MultiCaseResults
from HbO2.postprocess.tissuevolumes import TissueVolumeResults
from utilities.arguments import add_single_or_multi_case_group


def linregress_string(x, y):
    slope, intercept, rvalue, pvalue, stderr = linregress(x, y)
    return 'r = {:10.5g}, R^2 = {:10.5g}, p = {:7.5g}'.format(rvalue, rvalue**2, pvalue)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_single_or_multi_case_group(parser)
    args = parser.parse_args()
    multi_case = True if args.caseNames else False
    if multi_case:
        results = MultiCaseResults([TissueVolumeResults(case_path)
                                    for case_path in args.caseNames],
                                   combine_op=np.hstack)
    else:
        case_path = args.caseName
        results = TissueVolumeResults(case_path)
    print "Mean topol. radii:     ", np.mean(results.topological_tissue_radii())
    print "Mean func. radii:      ", np.mean(results.functional_tissue_radii())
    print "St. dev. topol. radii: ", np.std(results.topological_tissue_radii())
    print "St. dev. func. radii:  ", np.std(results.functional_tissue_radii())
    print
    print "Mean topol. volumes:     ", np.mean(results.topological_tissue_volumes())
    print "Mean func. volumes:      ", np.mean(results.functional_tissue_volumes())
    print "St. dev. topol. volumes: ", np.std(results.topological_tissue_volumes())
    print "St. dev. func. volumes:  ", np.std(results.functional_tissue_volumes())
    print "Sum topol. volumes:      ", np.sum(results.topological_tissue_volumes())
    print "Sum func. volumes:       ", np.sum(results.functional_tissue_volumes())
    print
    print "Functional vs topological radii: ", \
        linregress_string(results.functional_tissue_radii(), results.topological_tissue_radii())

    print "Functional radii vs mean hb:     ", \
        linregress_string(results.functional_tissue_radii(), results.hb_mean_on_edge())
    print "Functional radii vs hb drop:     ", \
        linregress_string(results.functional_tissue_radii(), results.hb_mean_drop())
    print "Functional radii vs RBC flow:    ", \
        linregress_string(results.functional_tissue_radii(), results.rbc_flow())
    print "Functional radii vs RBC vel.:    ", \
        linregress_string(results.functional_tissue_radii(), results.rbc_velocity())
    print "Functional radii vs hematocrit:  ", \
        linregress_string(results.functional_tissue_radii(), results.hematocrit())
    print "Functional radii vs transit time:", \
        linregress_string(results.functional_tissue_radii(), results.mean_arrival_transit_time())
    print "Functional radii vs path length: ", \
        linregress_string(results.functional_tissue_radii(), results.mean_arrival_path_length())
    print

    print "Topological radii vs mean hb:     ", \
        linregress_string(results.topological_tissue_radii(), results.hb_mean_on_edge())
    print "Topological radii vs hb drop:     ", \
        linregress_string(results.topological_tissue_radii(), results.hb_mean_drop())
    print "Topological radii vs RBC flow:    ", \
        linregress_string(results.topological_tissue_radii(), results.rbc_flow())
    print "Topological radii vs RBC vel.:    ", \
        linregress_string(results.topological_tissue_radii(), results.rbc_velocity())
    print "Topological radii vs hematocrit:  ", \
        linregress_string(results.topological_tissue_radii(), results.hematocrit())
    print "Topological radii vs transit time:", \
        linregress_string(results.topological_tissue_radii(), results.mean_arrival_transit_time())
    print "Topological radii vs path length: ", \
        linregress_string(results.topological_tissue_radii(), results.mean_arrival_path_length())
    print

