#!/usr/bin/env python
"""
Prepare an OpenFOAM case for an Lagrangian simulation.
Adapt the dictionaries for the mesh creation,
for the creation of RBC cellSets and modifies the dictionary
regionProperties.

Usage: prepareLagrangianCase.py
"""

import fileinput
import re
import shutil

from HbO2.setup.simulationParameters import IOHbO2SimulationParametersAxisymmetric
from preprocessing.changeEulerianMeshSize import changeMeshSize


blockMeshDict       = "domain/constant/polyMesh/blockMeshDict"
lagBlockMeshDict    = "lagrangian.org/constant/polyMesh/blockMeshDict"
sampleDict_x        = "domain/system/sampleDict.x_axis"
configVars          = "config_vars"
lagTopoSetDict      = "lagrangian.org/system/topoSetDict"

params = IOHbO2SimulationParametersAxisymmetric('.')
params.update_files()
nRBC = 1

# Change mesh
meshDict = changeMeshSize(
    params['domainLength'], [0,0,0],
    params['radiusTissueLeft'], params['radiusTissueRight'],
    params['radiusPlasma'], params['radiusWall'],
    params['dx'], params['dy'], params['nyTissue'],
    blockMeshDict)

# Create the set names
RBCRegionNames = ['RBC%i' % i for i in range(nRBC)]
lagMeshNames   = ['lagrangian%i' % i for i in range(nRBC)]

# Modify regionProperties and config_vars
regionPath = "domain/constant/regionProperties"
with open(regionPath, 'r') as f:
    for line in fileinput.input(regionPath, inplace=1):
        line = line.rstrip('\n')
        line = re.sub(r'^(\s*RBCRegionNames ).*',
                      '\\1( ' + ' '.join(RBCRegionNames) + ' );',
                      line)
        line = re.sub(r'^(\s*lagMeshNames ).*',
                      '\\1( ' + ' '.join(lagMeshNames) + ' );',
                      line)
        print line

with open(configVars, 'r') as f:
    for line in fileinput.input(configVars, inplace=1):
        line = line.rstrip('\n')
        line = re.sub(r'^(\s*RBCDirs).*',
                      '\\1=("' + '" "'.join(RBCRegionNames) + '");',
                      line)
        line = re.sub(r'^(\s*lagCaseDirs).*',
                      '\\1=("' + '" "'.join(lagMeshNames) + '");',
                      line)
        print line

print "Finished replacing entries in configuration files"

##############################################################
# Change Lagrangian case
##############################################################

lag_LDomain = 1.5*params['RBCLength']
# ensure that splitMeshRegions produces a mesh with the correct patches
# by forcing the radius of the lagrangian mesh to be larger than the RBC radius
# by (at least) one cell
lagMeshRadius = params['radiusRBC'] + 2*params['dyLag']
lagMeshDict = changeMeshSize(
    lag_LDomain, [-lag_LDomain/2, 0.0, 0.0],
    params['radiusTissueLeft'], params['radiusTissueRight'],
    lagMeshRadius, params['radiusWall'],
    params['dxLag'], params['dyLag'], 0,
    lagBlockMeshDict)

# Change STL file name in topoSetDict and define two 'outsidePoints', one at the lower
# corner, the other one at the upper right corner (this is safer if the Lagrangian box
# is only slightly bigger than the RBC).

if 'STLFile' in params:
    x_out1 = -0.49*lag_LDomain;
    y_out1 = 0.05*params['radiusPlasma'];
    z_out1 = 0.0;
    x_out2 = 0.49*lag_LDomain;
    y_out2 = 0.95*params['radiusPlasma'];
    z_out2 = 0.0;
    with open(lagTopoSetDict) as f:
        for line in fileinput.input(lagTopoSetDict, inplace=1):
            line = line.rstrip('\n')
            line = re.sub(r'^(\s*file\s*) "[^"]*";',
                          '\\1 "' + params['STLFile'] + '";',
                          line)
            line = re.sub(r'^(\s+outsidePoints\s+) [^;]+;',
                          '\\1 ((' + str(x_out1) + ' ' + str(y_out1) + ' ' + str(z_out1) + ') ' + \
                          '(' + str(x_out2) + ' ' + str(y_out2) + ' ' + str(z_out2) + '));',
                          line)

            print line

##############################################################
# Post-processing: adapt sampleDicts to mesh dimensions
##############################################################

# backup sampleDict_x
shutil.copyfile(sampleDict_x, sampleDict_x + '.bak')
# os.system("cp %s %s.bak" % (sampleDict_x, sampleDict_x))

# substitute domain half-length in sampleDict_x
with open(sampleDict_x) as f:
    for line in fileinput.input(sampleDict_x, inplace=1):
        line = line.rstrip('\n')
        # match expression in brackets after start
        m = re.search(r'^\s*(start|end)\s*\(([^)]+)\);', line)
        if (m != None):
            x_bound = params['domainLength']
            # reverse sign of x_bound if the entry starts with end
            if (m.group(1) == 'start'):
                x_bound = 0.0
            # split expression in brackets and replace the first entry
            s = [float(x) for x in m.group(2).split()]
            s[0] = float(x_bound)
            line = '        %s   ( % .7e % .7g % .7g );' % \
                   (m.group(1).ljust(5), s[0], s[1], s[2])

        print line

print "Finished preparing Lagrangian case."
