# Add and parse parser option for LD/RBC velocity study.

import argparse

def addOptions(parser):
    parser.add_argument('--LD', '-l', type=float, help='Linear density', default = 0.3)
    parser.add_argument('--vRBC', '-v', type=float, help='RBC velocity', default = 1e-3)
    parser.add_argument('-r', type=float, help='Radius for PO2 calculation', \
                        default = 17.6e-6)
    parser.add_argument('--cortex', '-c', help='Use cortex geometry (default)', \
                        default=True, action='store_true')
    parser.add_argument('--glomerulus', '-g', help='Use glomerulus geometry (precedes -c)', \
                        default=False, action='store_true')

def parseOptions(args):
    LD     = args.LD
    vRBC   = args.vRBC
    radius = args.r
    useGlomerulus = args.glomerulus
    useCortex = args.cortex

    if useGlomerulus:
        useCortex = False

    return {'LD': LD, 
            'vRBC': vRBC, \
            'radius': radius, \
            'useGlomerulus': useGlomerulus, \
            'useCortex': useCortex}


