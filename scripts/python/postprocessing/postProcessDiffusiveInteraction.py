#!/usr/bin/env python
"""
Postprocess simulations of diffusive interaction between capillaries
and compare results to experiments
"""

import argparse

from HbO2.COSH.integrate import makeDiffusiveInteractionIntegrators
from HbO2.postprocess.diffusiveinteraction import DiffusiveInteractionPostProcessor


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--nAverage', '-n', type=int, default=10,
                    help='Number of RBC paths used for averaging')
    args = parser.parse_args()
    integrators = makeDiffusiveInteractionIntegrators('.')
    pp = DiffusiveInteractionPostProcessor('.', integrators)
    pp.n_average = args.nAverage
    pp.write_results()
    print 'Wrote results to {:s}'.format(pp.output_file_name)
    print '{:-<{}}'.format('', 80)
    print '{:<{}s} {:.3e}'.format('Final x-coordinate:', 25, pp.finalXCoord())
    print '{:<{}s} {:g}'.format('Initial Hb mean:', 25, pp.initialHbMean())
    print '{:<{}s} {:g}'.format('Initial Hb difference:', 25, pp.initialHbDifference())
    print '{:-<{}}'.format('', 80)
    print '{:<{}} {:7.5g}'.format('Final Hb mean from simulation:', 40,
                                  pp.finalHbMeanFromSimul())
    print '{:<{}} {:7.5g}'.format('Final Hb mean from Krogh model:', 40,
                                  pp.finalHbMeanFromModel('krogh'))
    print '{:<{}} {:7.5g}'.format('Final Hb mean from simple model:', 40,
                                  pp.finalHbMeanFromModel('simple'))
    print '{:-<{}}'.format('', 80)
    print '{:<{}} {:7.5g}'.format('Absolute error in final Hb mean from model:', 40,
                                  pp.absModelErrorInFinalHbMean('krogh'))
    print '{:<{}} {:7.5%}'.format('Relative error in final Hb mean from model:', 40,
                                  pp.relModelErrorInFinalHbMean('krogh'))
    print '{:-<{}}'.format('', 80)
    print '{:<{}} {:7.5g}'.format('Final Hb difference from simulation:', 40,
                                  pp.finalHbDifferenceFromSimul())
    print '{:<{}} {:7.5g}'.format('Final Hb difference from Krogh model:', 40,
                                  pp.finalHbDifferenceFromModel('krogh'))
    print '{:<{}} {:7.5g}'.format('Final Hb difference from simple model:', 40,
                                  pp.finalHbDifferenceFromModel('simple'))
    print '{:<{}} {:7.5g}'.format('Final Hb difference from linearized model:', 40,
                                  pp.finalHbDifferenceFromModel('linearized'))
    print '{:-<{}}'.format('', 80)
    print '{:<{}} {:7.5g}'.format('Absolute error in final Hb difference from Krogh model:', 60,
                                  pp.absModelErrorInFinalHbDifference('krogh'))
    print '{:<{}} {:7.5g}'.format('Absolute error in final Hb difference from simple model:', 60,
                                  pp.absModelErrorInFinalHbDifference('simple'))
    print '{:<{}} {:7.5g}'.format('Absolute error in final Hb difference from linearized model:', 60,
                                  pp.absModelErrorInFinalHbDifference('linearized'))
    print '{:-<{}}'.format('', 80)
    print '{:<{}} {:7.5%}'.format('Relative error in final Hb difference from Krogh model:', 60,
                                  pp.relModelErrorInHbDifferenceDrop('krogh'))
    print '{:<{}} {:7.5%}'.format('Relative error in final Hb difference from simple model:', 60,
                                  pp.relModelErrorInHbDifferenceDrop('simple'))
    print '{:<{}} {:7.4%}'.format('Relative error in final Hb difference from linearized model:', 60,
                                  pp.relModelErrorInHbDifferenceDrop('linearized'))
