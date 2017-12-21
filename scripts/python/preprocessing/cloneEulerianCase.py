#!/usr/bin/env python

import os
import shutil
from optparse import OptionGroup
from os import path, system, chdir, mkdir

from PyFoam.Applications.PyFoamApplication import PyFoamApplication
from PyFoam.Error import error,warning
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory


class CloneEulerianCase(PyFoamApplication):
    def __init__(self,args=None):
        description="""\
Clone an Eulerian case. The directories 'domain' and
'lagrangian.org' are copied. The shell scripts to run the simulation
are also copied.\n
\n
Most files are copied using symlinks. The following\n
files are replaced by a copy of the original.\n
\n
domain/0.org\n
lagrangian.org/0.org\n
domain/system/controlDict\n
domain/system/fvSchemes\n
domain/system/fvSolution\n
domain/constant/polyMesh/blockMeshDict\n
domain/constant/polyMesh/blockMeshDict\n
config_vars\n
geometricData\n
initialConditions\n
simParams.json\n
.gitignore"""
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

        # copy domain and lagrangian.org
        mkdir(dName)
        chdir(dName)
        domainSrc=SolutionDirectory('../%s/%s' % (sName, 'domain'),
                              archive=None,
                              paraviewLink=False,
                              addLocalConfig=True,
                              parallel=self.opts.parallel)
        if self.parser.getOptions().symlink:
            domainSrc.symlinkCase(
                'domain',
                followSymlinks=True,
                maxLevel=1,
                relPath=True
            )
        else:
            domainSrc.cloneCase(
                'domain',
                followSymlinks=False
            )

        print 'Copied to %s/%s'  % (dName, 'domain.')

        lagSrc=SolutionDirectory('../%s/%s' % (sName, 'lagrangian.org'),
                              archive=None,
                              paraviewLink=False,
                              addLocalConfig=True,
                              parallel=self.opts.parallel)
        if self.parser.getOptions().symlink:
            lagSrc.symlinkCase(
                'lagrangian.org',
                followSymlinks=True,
                maxLevel=1,
                relPath=True
            )
        else:
            lagSrc.cloneCase(
                'lagrangian.org',
                followSymlinks=False
            )
        print 'Copied to %s/%s'  % (dName, 'lagrangian.org.')
        
        # replace domain.foam by <eulerianCaseName>.foam
        system('mv domain/domain.foam domain/%s.foam' % dName)

        # replace list of files by their original
        if self.parser.getOptions().symlink:
            replaceFiles = ['domain/system/controlDict',
                            'domain/system/fvSchemes',
                            'domain/system/fvSolution',
                            'domain/system/RBCInletDict',
                            'domain/constant/polyMesh/blockMeshDict']

            for f in replaceFiles:
                system('pyFoamSymlinkToFile.py %s' % f)

            system('pyFoamSymlinkToFile.py domain/0/*')

        # creates list of files to copy or to link to
        copyFiles = ['config_vars',
                     'geometricData',
                     'initialConditions',
                     'simParams.json',
                     '.gitignore',
                     'domain/.gitignore',
                     'lagrangian.org/.gitignore']

        symlinkFiles = []
        for f in os.listdir('../%s' % sName):
            if f.endswith('.stl'):
                symlinkFiles.append(f)

        optionDependentFiles = \
                    ['Allrun',
                     'AllrunRestart',
                     'cleanCase',
                     'createMesh',
                     'makeLinks',
                     'prepareEulerianCase',
                     'prepareLagrangianCases',
                     'removeLinks',
                     'runCase']

        if self.parser.getOptions().symlink:
            symlinkFiles.extend(optionDependentFiles)
        else:
            copyFiles.extend(optionDependentFiles)

        for f in copyFiles:
            srcName = os.path.join('..', sName,f)
            if os.path.isdir(srcName):
                if os.path.isdir(f):
                    shutil.rmtree(f)
                try:
                    shutil.copytree(srcName, f)
                except IOError:
                    print 'Could not copy the directory {}'.format(srcName)
            else:
                try:
                    shutil.copy(srcName, f)
                except IOError:
                    print 'Could not copy the directory {}'.format(srcName)

        for f in symlinkFiles:
            system('ln -s ../%s/%s %s' % (sName, f, f))


if __name__ == "__main__":
    CloneEulerianCase()


