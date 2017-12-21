#!/usr/bin/env python

import os

from PyFoam.Applications.PyFoamApplication import PyFoamApplication
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory

from HbO2.setup.case import copy_const_fields_first_to_last
from postprocessing.loadPO2Simulation import loadSimParams


class MyCopyLastToFirst(PyFoamApplication):
    def __init__(self,args=None):
        description="""\
Copies the contents of the last time-step of the source case to the
first time-step of the destination case (thus using it as initial
conditions). 

Additionally, copies the Lagrangian meshes and the RBC meshes of the 
source directory to the destination case."""
        PyFoamApplication.__init__(self,
                                   args=args,
                                   description=description,
                                   usage="%prog [options] <source caseDirectory> <destination caseDirectory>",
                                   interspersed=True,
                                   changeVersion=False,
                                   nr=2)

    def run(self):
        if len(self.parser.getArgs())!=2:
            error("Need two arguments.",len(self.parser.getArgs()),"found")

        # get directory names
        sName=self.parser.getArgs()[0]
        dName=self.parser.getArgs()[1]

        sDomain=sName + '/domain'
        dDomain=dName + '/domain'

        source=SolutionDirectory(sDomain,archive=None,paraviewLink=False)
        dest=SolutionDirectory(dDomain,archive=None,paraviewLink=False)

        sDir=source[-1]
        dDir=dest[0]

        # build the names of the Lagrangian and RBCs folders
        params = loadSimParams(sName)
        nRBC = params['nRBC']
        lagDirs = []
        RBCDirs = []
        for i in range(nRBC):
            lagDirs.append('lagrangian%i' % i)
            RBCDirs.append('RBC%i' % i)

        # in source, copy constant fields from the first to the last directory
        copy_const_fields_first_to_last(sDomain)

        # copy fields from the last directory of the source to the
        # first time step of the source
        copied=dDir.copy(sDir,
                         include='*',exclude=[],
                         overwrite=True,
                         mustExist=False,
                         purge=False)

        print 'Copied fields from %s/%s to %s/%s' % (sDomain, \
                sDir.baseName(), dDomain, dDir.baseName())

        # copy the folders 'lagrangian?'
        for f in lagDirs:
            sLag   = '%s/%s' % (sName,f)
            os.system('cp -r %s %s' % (sLag, dName))

        # copy 'points' files from the last time step of the source
        # domain to the 'constant/polyMesh' directory of the destination
        # and to the '0' directory.
        # Does this both for the Lagrangian and the RBC meshes.
        for (lagDir,RBCDir) in zip(lagDirs, RBCDirs):
            sLag = '%s/%s/%s' % (sDomain,sDir.baseName(),lagDir)
            dLag = '%s/%s/%s' % (dName, lagDir,'constant')
            d0Lag = '%s/%s/%s/' % (dDomain, '0', lagDir)
            os.system('cp %s/polyMesh/points* %s/polyMesh' % (sLag, dLag))
            os.system('cp %s/polyMesh/points* %s/polyMesh' % (sLag, d0Lag))

            sRBC = '%s/%s/%s' % (sDomain,sDir.baseName(),RBCDir)
            dRBC = '%s/%s/%s/%s' % (dName, lagDir,'constant', 'RBC')
            d0RBC = '%s/%s/%s' % (dDomain, '0', RBCDir)
            os.system('cp %s/polyMesh/points* %s/polyMesh' % (sRBC, dRBC))
            os.system('cp %s/polyMesh/points* %s/polyMesh' % (sRBC, d0RBC))

        print 'Copied Lagrangian meshes'

        # Copy fields defined on Lagrangian meshes and RBCs

        # Copy PO2_lag from last time step of the source to the
        # Lagrangian fields of the first time step of the destination
        for f in lagDirs:
            sLag = '%s/%s/%s' % (sDomain,sDir.baseName(),f)
            dLag = '%s/%s/%s' % (dName, f,'0')
            os.system('cp %s/PO2_lag* %s' % (sLag, dLag))
        print 'Copied PO2_lag'

        # Copy Hb from last time step of the source to the
        # RBC fields of the first time step of the destination
        for (lagDir,RBCDir) in zip(lagDirs, RBCDirs):
            sRBC = '%s/%s/%s' % (sDomain,sDir.baseName(),RBCDir)
            dRBC = '%s/%s/%s/%s' % (dName, lagDir,'0', 'RBC')
            os.system('cp %s/Hb* %s' % (sRBC, dRBC))
        print 'Copied Hb'


if __name__ == "__main__":
    MyCopyLastToFirst()


