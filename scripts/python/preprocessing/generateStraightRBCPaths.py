#!/usr/bin/env python
"""
Generate RBC paths in straight capillaries with a given linear density distribution
"""

import argparse

from HbO2.setup.rbcPaths import RBCPathGenerator


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    pathGenerator = RBCPathGenerator.fromparser(parser)
    pathGenerator.writeAll()
