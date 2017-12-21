#!/usr/bin/env python
"""
Replace patterns in files.

Patterns are read from a file. By default, all files in the current directory
are non-recursively traversed. All list of files can be given as an argument,
or the current directory can be recursively traversed.
"""

import argparse
import fileinput
import os
import re


def read_patterns(pattern_file):
    patterns = {}
    with open(pattern_file, 'r') as f:
        for line in f:
            l = line.split()
            patterns[l[0]] = l[1]
    return patterns


def read_filelist(filepath):
    filelist = []
    with open(filepath, 'r') as f:
        for line in f:
            filelist.append(line.rstrip('\n'))
    return filelist


def replace(filepath, patterns):
    for line in fileinput.input(filepath, inplace=1, backup=""):
        line = line.rstrip('\n')
        for pattern, after in patterns.iteritems():
            line = re.sub(pattern, after, line)
        print line

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', help='Recursive', action='store_true')
    parser.add_argument('--filelist', '-f', help='File name with file list')
    parser.add_argument('--patterns', '-p', help='File name with replace patterns')
    args = parser.parse_args()
    recursive = args.r
    filelist_path = args.filelist
    pattern_file = args.patterns

    patterns = read_patterns(pattern_file)

    if filelist_path:
        filelist = read_filelist(filelist_path)
        for f in filelist:
            replace(f, patterns)
            print 'Replaced pattern in file {}.'.format(f)
    elif recursive:
        for root, dirs, files in os.walk('.'):
            for f in files:
                if not f.endswith('.bak'):
                    filepath = os.path.join(root, f)
                    replace(filepath, patterns)
                    print 'Replaced pattern in file {}.'.format(filepath)
    else:
        for f in os.listdir('.'):
            if os.path.isfile(f) and not f.endswith('.bak'):
                replace(f, patterns)
                print 'Replaced pattern in file {}.'.format(f)

    print 'Finished replacing patterns.'
