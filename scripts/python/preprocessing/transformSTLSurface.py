#!/usr/bin/env python
#
# Translate a STL surface and writes the result to stdout.
#
# Usage: translateSTLSurface.py <STL_surface> <x> <y> <z>

import argparse

from HbO2.parse import STLUtils


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', \
                    help='path to surface')
parser.add_argument('--scale', '-s', type=float, help='scaling factor', \
                        default=1.0)
parser.add_argument('--tx', '-x', type=float, \
                        help='x-component of translation vector', \
                        default=0.0)
parser.add_argument('--ty', '-y', type=float, \
                        help='y-component of translation vector', \
                        default=0.0)
parser.add_argument('--tz', '-z', type=float, \
                        help='z-component of translation vector', \
                        default=0.0)
args = parser.parse_args()
filePath = args.file
scaling_factor = args.scale
tx       = args.tx
ty       = args.ty
tz       = args.tz

STLUtils.transformSurface(filePath, scaling_factor, [tx, ty, tz])
