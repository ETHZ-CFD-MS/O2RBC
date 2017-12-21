#!/usr/bin/env python

import os
import shutil
from optparse import OptionGroup
from os import path, system, chdir, mkdir

from PyFoam.Applications.PyFoamApplication import PyFoamApplication
from PyFoam.Error import error,warning
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory


class CloneGraphCase(PyFoamApplication):
    def __init__(self,args=None):
        description="""
Clone an graph case in the CBF framework. The directories 'domain'
'RBC' and 'RBCSource' are copied. The shell scripts to run the simulation
are also copied.\n
Files can be copied using symlinks.
"""
        PyFoamApplication.__init__(self,
                                   args=args,
                                   description=description,
                                   usage="%prog [options] <source caseDirectory> <destination caseDirectory>",
                                   interspersed=True,
                                   changeVersion=False,
                                   nr=2)

    def addOptions(self):
        behave=OptionGroup(self.parser,
                           "Behaviour")
        self.parser.add_option_group(behave)
        behave.add_option("--symlink",
                          action="store_true",
                          dest="symlink",
                          default=False,
                          help="Clone using symlink")
        behave.add_option("--parallel",
                          action="store_true",
                          dest="parallel",
                          default=False,
                          help="Clone the processor-directories")
        behave.add_option("--force",
                          action="store_true",
                          dest="force",
                          default=False,
                          help="Overwrite destination if it exists")
    def run(self):
        if len(self.parser.getArgs())!=2:
            error("Need two arguments.",len(self.parser.getArgs()),"found")

        # get directory names
        sName=self.parser.getArgs()[0]
        dName=self.parser.getArgs()[1]

        # check if the destination directory exists
        if path.exists(dName):
            if self.parser.getOptions().force:
                warning("Replacing",dName,"(--force option)")
            elif path.exists(path.join(dName,"system","controlDict")):
                error("Destination",dName,"already existing and a Foam-Case")
            elif path.isdir(dName):
                dName=path.join(dName,path.basename(sName))
                if path.exists(dName) and not self.parser.getOptions().force:
                    error(dName,"already existing")
        elif not path.exists(path.dirname(dName)):
            warning("Directory",path.dirname(dName),"does not exist. Creating")

        mkdir(dName)
        chdir(dName)

        self.cloneCase(os.path.join('..', sName, 'domain'), 'domain')
        self.cloneCase(os.path.join('..', sName, 'RBC'), 'RBC')

        RBCSourcePath = os.path.join('..', sName, 'RBCSource')
        if os.path.isdir(RBCSourcePath):
            self.cloneCase(RBCSourcePath, 'RBCSource')

        sampleRBCPath = os.path.join('..', sName, 'sampleRBC')
        if os.path.isdir(sampleRBCPath):
            self.cloneCase(sampleRBCPath, 'sampleRBC')

        # replace domain.foam by <eulerianCaseName>.foam
        shutil.move(os.path.join('domain', 'domain.foam'),
                    os.path.join('domain', '{}.foam'.format(dName)))

        # replace list of files by their original
        if self.parser.getOptions().symlink:
            replaceFiles = ['domain/system/controlDict',
                            'domain/system/fvSchemes',
                            'domain/system/fvSolution',
                            'domain/constant/polyMesh/blockMeshDict']

            for f in replaceFiles:
                system('pyFoamSymlinkToFile.py %s' % f)

            system('pyFoamSymlinkToFile.py domain/0.org/*')

        # creates list of files to copy or to link to
        copyFiles = ['config_vars',
                     'initialConditions',
                     'geometricData']

        symlinkFiles = []
        for f in os.listdir(os.path.join('..', sName)):
            if f.endswith('.stl'):
                symlinkFiles.append(f)

        optionDependentFiles = \
                    ['Allrun',
                     'AllrunRestart',
                     'cleanCase',
                     'prepareEulerianCase',
                     'createRBCMesh',
                     'makeLinks',
                     'removeLinks',
                     'restartSimulation',
                     'runCase',
                     'graphDict.pkl',
                     'RBCPathParams.json',
                     'simParams.json']

        if self.parser.getOptions().symlink:
            symlinkFiles.extend(optionDependentFiles)
        else:
            copyFiles.extend(optionDependentFiles)

        for f in copyFiles:
            copyIfFileExists(os.path.join('..', sName, f), f)

        for f in symlinkFiles:
            os.symlink(os.path.join('..', sName, f), f)

    def cloneCase(self, srcPath, dstName):
        dName=self.parser.getArgs()[1]
        srcCase=SolutionDirectory(srcPath,
                                  archive=None,
                                  paraviewLink=False,
                                  addLocalConfig=True,
                                  parallel=self.opts.parallel)
        if self.parser.getOptions().symlink:
            srcCase.symlinkCase(
                dstName,
                followSymlinks=True,
                maxLevel=1,
                relPath=True
            )
        else:
            srcCase.cloneCase(
                dstName,
                followSymlinks=False
            )
        print 'Copied to {:s}'.format(os.path.join(dName, dstName))


def copyIfFileExists(src, dst):
    """
    Copy the file with path src to dst if the file exists.

    Args:
        src (str): source file
        dst (str): destination file (whole path)
    """
    if os.path.exists(src) and os.path.isfile(src):
        shutil.copy(src, dst)

if __name__ == "__main__":
    CloneGraphCase()



