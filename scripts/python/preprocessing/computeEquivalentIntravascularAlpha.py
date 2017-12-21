#!/usr/bin/env python
"""
Determines the constant solubility coefficient in the RBC, plasma and wall so that
the resulting analytical intravascular coefficient is equal to that with the given
coefficients.
"""

import numpy as np

from HbO2.model.coefficients import intravascularResistanceAnalytical
from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric


def equivalent_alpha(simParams, ivr):
    """
    Compute the equivalent alpha

    Args:
        simParams (IOHbO2ParametersAxisymmetric): simulation parameters

    Returns:
        (float) alpha
    """
    LD = simParams['LDMean']
    kappaO2RBC = simParams['kappaO2RBC']
    kappaO2Plasma = simParams['kappaO2Plasma']
    kappaO2Wall = simParams['kappaO2Wall']
    Rrbc = simParams['radiusRBC']
    Rp = simParams['radiusPlasma']
    Rw = simParams['radiusWall']
    return 1./(2*np.pi*LD*ivr)*(1./(4*kappaO2RBC)
                              + np.log(Rp/Rrbc)/kappaO2Plasma
                              + np.log(Rw/Rp)/kappaO2Wall)


if __name__ == '__main__':
    simParams = IOHbO2ParametersAxisymmetric('.')
    ivr_target = intravascularResistanceAnalytical(simParams)
    print 'Equivalent alpha = {:g}'.format(equivalent_alpha(simParams, ivr_target))
