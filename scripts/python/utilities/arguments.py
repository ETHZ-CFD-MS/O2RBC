"""
Utility functions for argument parsing
"""

import argparse


def add_case_argument(parser):
    """
    Add argument for a case name a parser

    Args:
        parser (argparse.ArgumentParser): argument parser
    """
    parser.add_argument('--caseName', '-c', type=str, help='Case name', default='.')


def add_multi_case_argument(parser):
    """
    Add required arguments for minimal and maximal time to parser

    Args:
        parser (argparse.ArgumentParser): argument parser
    """
    parser.add_argument('--caseNames', type=str, nargs='+', help='Case names')


def add_single_or_multi_case_group(parser):
    """
    Add a mutually exclusive group for one or multiple cases.

    Args:
        parser (argparse.ArgumentParser): argument parser
    """
    case_group = parser.add_mutually_exclusive_group()
    add_case_argument(case_group)
    add_multi_case_argument(case_group)


def add_min_max_time_argument(parser):
    """
    Add required arguments for minimal and maximal time to parser

    Args:
        parser (argparse.ArgumentParser): argument parser
    """
    parser.add_argument('--minTime', type=float, help='Minimal time', default=-1e300)
    parser.add_argument('--maxTime', type=float, help='Maximal time', default=1e300)
