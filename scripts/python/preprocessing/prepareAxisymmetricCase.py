#!/usr/bin/env python
"""
Prepare an OpenFOAM case for a simulation in an axisymmetric domain.
Adapt the dictionaries for the mesh creation.
"""

import fileinput
import re

from HbO2.setup.simulationParameters import IOHbO2SimulationParametersAxisymmetric
from preprocessing.changeEulerianMeshSize import changeMeshSize


domain_block_mesh_dict_path = "domain/constant/polyMesh/blockMeshDict"
rbc_block_mesh_dict_path = "sampleRBC/constant/polyMesh/blockMeshDict"
sample_dict_x_path = "domain/system/singleGraph.x_axis"

params = IOHbO2SimulationParametersAxisymmetric('.')
params.update_files()

# Change Eulerian mesh
changeMeshSize(
    params['domainLength'], [0, 0, 0],
    params['radiusTissueLeft'], params['radiusTissueRight'],
    params['radiusPlasma'], params['radiusWall'],
    params['dx'], params['dy'], params['nyTissue'],
    domain_block_mesh_dict_path)

# Change mesh for RBC
changeMeshSize(
    params['RBCLength'], [-params['RBCLength']/2, 0.0, 0.0],
    params['radiusTissueLeft'], params['radiusTissueRight'],
    params['radiusRBC'], params['radiusWall'],
    params['dxRBC'], params['dyRBC'], 0,
    rbc_block_mesh_dict_path)

# Post-processing: adapt sample dictionary to mesh dimensions
with open(sample_dict_x_path) as f:
    for line in fileinput.input(sample_dict_x_path, inplace=1, backup=""):
        line = line.rstrip('\n')
        # match expression in brackets after start
        m = re.search(r'^\s*(start|end)\s*\(([^)]+)\);', line)
        if m is not None:
            x_bound = params['domainLength']
            # reverse sign of x_bound if the entry starts with end
            if m.group(1) == 'start':
                x_bound = 0.0
            # split expression in brackets and replace the first entry
            s = [float(x) for x in m.group(2).split()]
            s[0] = float(x_bound)
            line = '        %s   ( % .7e % .7g % .7g );' % \
                   (m.group(1).ljust(5), s[0], s[1], s[2])
        print line
