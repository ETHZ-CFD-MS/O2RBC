#!/usr/bin/env python
#
# Plot an OpenFOAM mesh
#

import argparse
import matplotlib.pyplot as plt

from parse import readMesh
from plot.plotMesh import plotMesh
from plot.figureoptions import FigureOptions


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    figOptions = FigureOptions(parser)

    args            = parser.parse_args()
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    RBCStyle    = {'color' : '0.6',
                   'linewidth' :  0.6}
    domainStyle = {'color' : 'k',
                   'linewidth' : 0.2}

    # read the RBC mesh
    (pointsRBC, facesRBC) = readMesh.readMesh('lagrangian0', timePath=None,
            regionName='RBC')

    # read the Eulerian mesh
    (points, faces) = readMesh.readMesh('domain')

    # plot the Eulerian mesh on the top left
    # translate the RBC mesh
    xt = 5e-5
    points[:,0] = points[:,0] + xt
    plotMesh(points, faces, style=RBCStyle)

    # plot the domain mesh
    plotMesh(points, faces, style=domainStyle)

    # plt.gca().set_xlim([40, 60])
    figOptions.adjustAxes()
    # plt.axis('equal')
    # ax_pt = plt.gca().get_position().get_points()
    # print ax_pt
    # ratio = (ax_pt[1][1]-ax_pt[0][1])/(ax_pt[1][0]-ax_pt[0][0])
    xw = 8.5
    xc = 50
    plt.xlim([xc-xw/2, xc+xw/2])
    plt.ylim([0, 0.85*6])
    # plt.ylim([0, ratio*xw])
    plt.axis('off')

    plotName = 'plotOverlappingMeshes'
    # plotName = 'plotDomainMesh'
    # plotName = 'plotRBCMesh'
    # plt.savefig('%s.svg' % plotName, transparent=True)
    figOptions.saveFig(plotName)



