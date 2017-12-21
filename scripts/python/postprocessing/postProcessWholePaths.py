#!/usr/bin/env python
"""
Postprocess tissue volumes obtained from simulations in a graph
"""
from __future__ import division
import argparse
import numpy as np
from scipy.stats.stats import linregress

from HbO2.postprocess.factory.case import make_post_processor
from utilities.arguments import add_case_argument


def linregress_string(x, y):
    slope, intercept, rvalue, pvalue, stderr = linregress(x, y)
    return 'R^2 = {:g}, p = {:g}'.format(rvalue**2, pvalue)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_case_argument(parser)
    parser.add_argument('--settingsFile', default='postprocess.json')
    args = parser.parse_args()
    case_path = args.caseName
    settings_file = args.settingsFile
    threshold_low_hb = 0.2
    pp = make_post_processor(case_path, settings_file)

    sim_distal_hb = pp.distal_hb()
    sim_distal_mean = np.mean(pp.distal_hb())
    sim_distal_std = np.std(pp.distal_hb())
    sim_proximal_hb = pp.proximal_hb()
    sim_proximal_mean = np.mean(pp.proximal_hb())
    sim_proximal_std = np.std(pp.proximal_hb())
    sim_fraction_low = len(np.where(sim_distal_hb <= threshold_low_hb)[0])/len(sim_distal_hb)

    pp.integrator.use_topological_radii = True
    pp.integrator.compute()
    int_topol_distal_hb, topol_weights = pp.integrator.distal_hb_distribution()
    int_topol_distal_mean = pp.integrator.average_distal_hb()
    int_topol_distal_std = pp.integrator.std_distal_hb()
    int_topol_distal_fraction_low = pp.integrator.distal_hb_fraction_below(threshold_low_hb)

    pp.integrator.use_topological_radii = False
    pp.integrator.compute()
    int_func_distal_hb, func_weights = pp.integrator.distal_hb_distribution()
    int_func_distal_mean = pp.integrator.average_distal_hb()
    int_func_distal_std = pp.integrator.std_distal_hb()
    int_func_distal_fraction_low = pp.integrator.distal_hb_fraction_below(threshold_low_hb)

    # statistics for proximal and distal hemoglobin distribution
    print "{:<30s} {:g}".format("Mean simulated proximal hb:", sim_proximal_mean)
    print "{:<30s} {:g}".format("Std simulated proximal hb:", sim_proximal_std)
    print
    print "{:<30s} {:g}".format("Mean simulated distal hb:", sim_distal_mean)
    print "{:<30s} {:g}".format("Std simulated distal hb:", sim_distal_std)
    print
    print "{:<30s} {:g}".format("Mean int. distal hb (func):", int_func_distal_mean)
    print "{:<30s} {:g}".format("Std int. distal hb (func):", int_func_distal_std)
    print "{:<30s} {:g}".format("Mean int. distal hb (topol):", int_topol_distal_mean)
    print "{:<30s} {:g}".format("Std int. distal hb (topol):", int_topol_distal_std)
    print
    print "{:<30s} {:g}".format("Mean hb drop:", pp.mean_hb_drop())
    print "{:<30s} {:g}".format("Std hb drop:", pp.std_hb_drop())
    print "{:<30s} {:g}".format("CoV hb drop:", pp.coeff_variation_hb_drop())
    print
    print "{:<30s} {:g}".format("Sim. frac. S <= {:g}:".format(threshold_low_hb),
                                sim_fraction_low)
    print "{:<30s} {:g}".format("Int. frac. S <= {:g} (func):".format(threshold_low_hb),
                                int_func_distal_fraction_low)
    print "{:<30s} {:g}".format("Int. frac. S <= {:g} (topol):".format(threshold_low_hb),
                                int_topol_distal_fraction_low)
    print

    # statistics for transit times
    print "{:<30s} {:g}".format("Mean transit time:", pp.mean_transit_time())
    print "{:<30s} {:g}".format("Std transit time:", pp.std_transit_time())
    print "{:<30s} {:g}".format("CoV transit time:", pp.coeff_variation_transit_time())
    print

    # statistics for transit paths
    print "{:<30s} {:g}".format("Mean transit path length:", pp.mean_transit_path_length())
    print "{:<30s} {:g}".format("Std transit path length:", pp.std_transit_path_length())
    print "{:<30s} {:g}".format("CoV transit path length:", pp.coeff_variation_transit_path_length())
    print

    # correlations
    print "{:<40s} {:g}, p = {:e}".format("Sim. hb drop vs transit times:",
                                          *pp.corr_sim_hb_drop_vs_transit_time())
    print "{:<40s} {:g}, p = {:e}".format("Sim. hb drop vs transit paths:",
                                          *pp.corr_sim_hb_drop_vs_transit_path_length())
    print "{:<40s} {:g}, p = {:e}".format("Hb drop vs proximal hb:",
                                          *pp.corr_hb_drop_vs_proximal_hb())
    print "{:<40s} {:g}, p = {:e}".format("Hb slope transit time vs proximal hb:",
                                          *pp.corr_mean_hb_slope_transit_time_vs_proximal_hb())
    print "{:<40s} {:g}, p = {:e}".format("Hb slope path length vs proximal hb:",
                                          *pp.corr_mean_hb_slope_path_length_vs_proximal_hb())

    # relative reduction of COSH
    print
    print "{:<52s} {:8.3%}".format("COSH reduction by diff. inter.:",
                                1 - sim_distal_std/int_topol_distal_std)
    print "{:<52s} {:8.3%}".format("COSH reduction by cap. diff. inter.:",
                                1 - int_func_distal_std/int_topol_distal_std)
    print "{:<52s} {:8.3%}".format("Fraction of COSH reduction due to cap. diff. inter.:",
                                (int_topol_distal_std - int_func_distal_std)/
                                (int_topol_distal_std - sim_distal_std))
    print "{:<52s} {:8.3%}".format("Fraction of COSH reduction due to RBC diff. inter.:",
                                (int_func_distal_std - sim_distal_std)/
                                (int_topol_distal_std - sim_distal_std))
