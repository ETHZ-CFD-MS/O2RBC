#!/usr/bin/env python
#
# Copy constant fields required for the simulations from the time directory 0
# to the last time directory.
#
# Usage: copy_const_fields_first_to_last.py --help

import argparse
import os

from HbO2.setup.case import copy_const_fields_first_to_last

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('case', nargs='?', help='case directory', default='domain')
    parser.add_argument('-p', '--parallel', action='store_true')
    args = parser.parse_args()
    sName = args.case
    parallel = args.parallel

    if parallel:
        copied = []
        lastTime = ''
        for d in os.listdir(sName):
            if d.startswith('processor'):
                sNameProc = os.path.join(sName, d)
                copied, lastTime = copy_const_fields_first_to_last(sNameProc)
    else:
        copied, lastTime = copy_const_fields_first_to_last(sName)

    filename = os.path.basename(__file__)
    if copied:
        print '{:s}: Copied fields {:s} to time directory {:s}.'.format(filename,
                                                                        ', '.join(copied), lastTime)
    else:
        print '{:s}: No fields were copied.'.format(filename)
