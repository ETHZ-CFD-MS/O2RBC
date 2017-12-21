#!/usr/bin/env python
#
# Compute Krogh solution for PO2 in tissue, given PO2 of the
# capillary.
#

import numpy as np


kappa = 2.41e-9
alpha = 38.9 

def computeKroghSolution(r, P_cap, R_c, R_t, M):

    P = P_cap + M/(4*kappa*alpha) * (r**2 + R_c**2) \
              - M*R_t**2/(2*kappa*alpha) * np.log(r/R_c)
    print P

if __name__ == "__main__":
    
    n = 100 # number of grid points

    R_c   = 2e-6
    R_t   = 26e-6
    M     = 1830
    factor = 1.5

    P_cap = 37.8

    r = np.linspace(R_c, R_t, n)

    computeKroghSolution(R_t, P_cap, R_c, R_t, factor*M)


