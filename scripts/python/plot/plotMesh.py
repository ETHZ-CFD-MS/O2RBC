#!/usr/bin/env python
#
# Plot an OpenFOAM mesh
#

import argparse
import matplotlib.pyplot as plt

from parse import readMesh
from plot.figureoptions import FigureOptions


def plotMesh(points, faces, style=None):

    if style is None:
        style = {'edgecolor' : 'k',
                 'facecolor' : '0.8'}
    # convert to micrometers
    points = 1e6*points

    xPoints = points[:,0]
    yPoints = points[:,1]

    # plot vertices
    # plt.plot(xPoints, yPoints, 'k.')

    # plot faces
    xFaces = []
    yFaces = []
    for face in faces:
        xFaces.append([xPoints[i] for i in face])
        yFaces.append([yPoints[i] for i in face])

    # if style['facecolor'] is not None:
        # for xFace, yFace in zip(xFaces, yFaces):
            # plt.fill(xFace, yFace, **style)
    # else:
        # Faster code taken from
        # http://exnumerus.blogspot.ch/2011/02/how-to-quickly-plot-polygons-in.html
        # Seems not to work for the facecolor
    xList = []
    yList = []
    for xFace, yFace in zip(xFaces, yFaces):
        xList.extend(xFace)
        xList.append(None)
        yList.extend(yFace)
        yList.append(None)

    plt.plot(xList, yList, **style)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    figOptions = FigureOptions(parser)

    args            = parser.parse_args()
    figOptions.parseOptions(args)
    figOptions.applyOptions()

    # RBCStyle    = {'edgecolor' : 'k',
                   # 'facecolor' : None,
                   # 'linewidth' :  0.6,
                   # 'alpha'     :  0.5}
    # domainStyle = {'edgecolor' : 'k',
                   # 'facecolor' : None,
                   # 'linewidth' : 0.2}
    RBCStyle    = {'color' : '0.6',
                   'linewidth' :  0.6}
    domainStyle = {'color' : 'k',
                   'linewidth' : 0.2}

    (points, faces) = readMesh.readMesh('lagrangian0', timePath=None,
            regionName='RBC')
    # translate the RBC mesh
    xt = 5e-5
    points[:,0] = points[:,0] + xt
    plotMesh(points, faces, style=RBCStyle)

    # plot the domain mesh
    (points, faces) = readMesh.readMesh('domain')
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


