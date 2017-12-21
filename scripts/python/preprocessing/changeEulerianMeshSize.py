#!/usr/bin/env python
"""
Change the dimensions of an axisymmetric mesh of an OpenFOAM case 
that is generated using a blockMeshDict file.

Usage: changeEulerianMeshSize.py --help
"""

import argparse
import re
import fileinput
import math
import numpy as np


def changeMeshSize(L_domain, origin, Rt_left, Rt_right, Rp, Rw, 
                   dx, dy, ny_t, blockMeshDict):
    """Changes a blockMeshDict file by substituting the given values.

    Args:
        L_domain: domain length (x-direction)
        origin:   coordinates of lower left corner
        Rt_left:  tissue radius on the left
        Rt_right: tissue radius on the right
        Rp:       plasma radius
        Rw:       wall radius
        dx:       mesh width in x direction
        dy:       mesh width in y-direction (in plasma and wall)
        ny_t:     number of grid cells in y-direction in tissue
        blockMeshDict: path to blockMeshDict file
    """

    # parameters
    alpha_deg = 5 # wedge angle in degrees

    # computed values
    alpha = alpha_deg*math.pi/180

    tissueVolumeGoal = volumeTissueAnalytical(Rt_left, Rt_right, Rw, L_domain, alpha)
    tissueMeshVolume = volumeTissueMesh(Rt_left, Rt_right, Rw, L_domain, alpha)
    delta = computeDeltaCompensation(Rt_left, Rt_right, Rw, L_domain, alpha)
    newTissueMeshVolume = volumeTissueMesh(Rt_left+delta, Rt_right+delta, Rw, L_domain, alpha)
    print "Adapted mesh tissue volume from {:g} to {:g}.".format(tissueMeshVolume, newTissueMeshVolume)

    x_l  = origin[0]
    x_r  = origin[0] + L_domain
    y_p  = origin[1] + Rp              *math.cos(alpha/2)
    y_w  = origin[1] + Rw              *math.cos(alpha/2)
    y_tl = origin[1] + (Rt_left +delta)*math.cos(alpha/2)
    y_tr = origin[1] + (Rt_right+delta)*math.cos(alpha/2)
    z_p  = origin[2] + Rp              *math.sin(alpha/2)
    z_w  = origin[2] + Rw              *math.sin(alpha/2)
    z_tl = origin[2] + (Rt_left +delta)*math.sin(alpha/2)
    z_tr = origin[2] + (Rt_right+delta)*math.sin(alpha/2)
    nx   = int(round(L_domain/dx))
    ny_p = int(round(Rp/dy))
    ny_w = int(round((Rw-Rp)/dy))

    # replace entries in blockMeshDict
    with open(blockMeshDict, 'r') as f:
        for line in fileinput.input(blockMeshDict, inplace=1):
            line = line.rstrip('\n')
            line = re.sub(r'^(x_l\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), x_l),
                          line)
            line = re.sub(r'^(x_r\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), x_r),
                          line)
            line = re.sub(r'^(y_p\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), y_p),
                          line)
            line = re.sub(r'^(y_w\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), y_w),
                          line)
            line = re.sub(r'^(y_tl\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), y_tl),
                          line)
            line = re.sub(r'^(y_tr\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), y_tr),
                          line)
            line = re.sub(r'^(z_p\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), z_p),
                          line)
            line = re.sub(r'^(z_w\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), z_w),
                          line)
            line = re.sub(r'^(z_tl\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), z_tl),
                          line)
            line = re.sub(r'^(z_tr\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), z_tr),
                          line)
            line = re.sub(r'^(mz_p\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), -z_p),
                          line)
            line = re.sub(r'^(mz_w\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), -z_w),
                          line)
            line = re.sub(r'^(mz_tl\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), -z_tl),
                          line)
            line = re.sub(r'^(mz_tr\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), -z_tr),
                          line)
            line = re.sub(r'^(nx\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), nx),
                          line)
            line = re.sub(r'^(ny_p\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), ny_p),
                          line)
            line = re.sub(r'^(ny_w\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), ny_w),
                          line)
            line = re.sub(r'^(ny_t\s+)[^;]*;',
                          lambda m: "%s%.9g;" % (m.group(1), ny_t),
                          line)
            print line

    print "Changed entries in file %s." % blockMeshDict

    return {'L_domain': L_domain, 
            'Rt_left' : Rt_left, 
            'Rt_right': Rt_right,
            'Rp'      : Rp, 
            'Rw'      : Rw, 
            'dx'      : dx, 
            'dy'      : dy, 
            'nx'      : nx,
            'ny_p'    : ny_p, 
            'ny_w'    : ny_w, 
            'ny_t'    : ny_t}

def replaceVariable(line, key, value):
    line = re.sub(r'^(x_l\s*)[^;]*;',
                  lambda m: '%s %f;' % (m.group(1), value),
                  line)
            
def volumeTissueMesh(Rt_left, Rt_right, Rw, L, alpha):
    return L*math.sin(alpha/2)*math.cos(alpha/2)* \
            (1./3.*(Rt_left**2 + Rt_left*Rt_right + Rt_right**2) - Rw**2)

def volumeTissueAnalytical(Rt_left, Rt_right, Rw, L, alpha):
    return L*(alpha/2)* \
            (1./3.*(Rt_left**2 + Rt_left*Rt_right + Rt_right**2) - Rw**2)

def computeDeltaCompensation(Rt_left, Rt_right, Rw, L, alpha):
    roots = np.roots([1,
                      Rt_left + Rt_right,
                     -(1./3.*(Rt_left**2 + Rt_left*Rt_right + Rt_right**2) - Rw**2) \
                     *(alpha/2/(math.sin(alpha/2)*math.cos(alpha/2)) - 1)])
    if (max(roots) < 0):
        print "No positive root for delta, returning 0!"
        return 0
    return max(roots)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--length', '-L', type=float, help='domain length')
    parser.add_argument('--Rt', '-R', type=float, help='tissue radius',
                        default=18e-6)
    parser.add_argument('--Rc', type=float, help='capillary radius',
                        default=3e-6)
    parser.add_argument('--dx', type=float, help='dx', default=2e-7)
    parser.add_argument('--dy', type=float, help='dy', default=2e-7)
    parser.add_argument('--ny_t', type=float,
                        help='number of cells in radial direction in the tissue',
                        default=100)
    parser.add_argument('-f', '--file',
                        help='path to blockMeshDict',
                        default='constant/polyMesh/blockMeshDict')

    args = parser.parse_args()

    L_domain = args.length
    Rt       = args.Rt
    Rc       = args.Rc
    dx       = args.dx
    dy       = args.dy
    ny_t     = args.ny_t
    blockMeshDict = args.file

    changeMeshSize(L_domain, Rt, Rc, dx, dy, ny_t, blockMeshDict)
