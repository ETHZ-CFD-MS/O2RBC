#!/usr/bin/env python
#
# Sample using OpenFOAM for all singleGraph files.
#

import argparse
import glob
import os
import shutil

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--overwrite', '-o', action='store_true')
args = parser.parse_args()
overwrite = args.overwrite

setFolderPath = 'postProcessing'

if overwrite:
    for dir in glob.glob(os.path.join(setFolderPath, 'sets*')):
        shutil.rmtree(dir)
    print 'Cleaned up old set directories'

# get all singleGraph files
sampleDicts = []
for f in os.listdir('domain/system'):
    if f == 'singleGraph' or f.startswith('singleGraph.'):
        sampleDicts.append(f)

has_skipped_folder = False
# for each sampleDict file, sample and move folders
for sampleDict in sampleDicts:
    # get the extension of the dictionary name
    extension = os.path.splitext(sampleDict)[1]
    if extension == '':
        extension = 'default'
    else:
        # remove dot
        extension = extension[1:]

    # create folder names
    oldSetName = 'sets'
    oldSetPath = os.path.join(setFolderPath, oldSetName)
    newSetName = oldSetName + '_' + extension
    newSetPath = os.path.join(setFolderPath, newSetName)

    if not os.path.exists(newSetPath):
        # sample using the OpenFOAM utility
        os.system('postProcess -case domain -dict {:s}'.format(os.path.join('system', sampleDict)))
        # move the sets directory
        if os.path.exists(os.path.join(setFolderPath, oldSetName)):
            shutil.move(oldSetPath, newSetPath)
    else:
        has_skipped_folder = True
        print 'Folder {:s} already present, skipping it.'.format(newSetName)

if has_skipped_folder:
    print 'The option -o/--overwrite can be used to overwrite old results'

