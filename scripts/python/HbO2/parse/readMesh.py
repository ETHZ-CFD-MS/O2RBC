#!/usr/bin/env python
"""
Functions for reading a mesh produced by OpenFOAM
"""

import numpy as np


def readMesh(caseName, timePath=None, regionName=None):
    """
    Reads an OpenFOAM mesh

    Args:
        caseName: Path to the OpenFOAM case
        timePath: Name of time folder, if the points of interest are in there
        regionName: Name of the mesh region, if any.
    Returns:
        tuple with points and faces
    """
    polyMesh = pathToPolyMesh(caseName, timePath=timePath, regionName=regionName)
    points = readPoints(polyMesh)
    faces  = readFaces(polyMesh)
    return (points, faces)


def readPoints(polyMesh):
    """
    Read the vertex coordinates in an OpenFOAM mesh

    Args:
        polyMesh: Path to the polyMesh folder

    Returns:
        numpy array with point coordinates
    """
    pathToPoints = '%s/points' % polyMesh
    # find the number of points
    (nPoints, lineNPoints) = readNumberItemsOFDict(pathToPoints)
    # parse the point positions
    points = np.zeros( (nPoints, 3) )
    with open(pathToPoints, 'r') as f:
        i = 0
        for line in f:
            i = i+1
            line = line.rstrip('\n')
            iPoint = i - (lineNPoints+2)
            if iPoint >= 0 and iPoint < nPoints:
                # remove brackets from line
                s = line.translate(None, '()')
                # copy coordinates into array
                points[iPoint,:] = [float(x) for x in s.split()]
    return points


def readFaces(polyMesh):
    """
    Read the face index of an OpenFOAM mesh

    Args:
        polyMesh Path to the polyMesh folder

    Returns:
        List of lists with points indices
    """

    pathToFaces = '%s/faces' % polyMesh
    # find the number of faces
    (nFaces, lineNFaces) = readNumberItemsOFDict(pathToFaces)
    # parse the faces
    faces = []
    with open(pathToFaces, 'r') as f:
        i = 0
        for line in f:
            i = i+1
            line = line.rstrip('\n')
            iFace = i - (lineNFaces+2)
            if iFace >= 0 and iFace < nFaces:
                # extract the number of vertices
                nv = line[0] # FIXME, not robust
                # remove brackets from line
                s = line[1:].translate(None, '()')
                # copy coordinates into array
                faces.append([int(x) for x in s.split()])
    return faces

def readNumberItemsOFDict(pathToFile):
    """
    Find the number of items in an OpenFOAM dictionary

    Args:
        pathToFile: Path to dictionary

    Returns:
        Tuple with the number of items and the line number where it was found
    """
    nItems = -1
    with open(pathToFile, 'r') as f:
        i = 0
        for line in f:
            i = i+1
            # if the line is composed of digits only, it is the number
            # of items
            line = line.rstrip('\n')
            if line.isdigit():
                nItems = float(line)
                lineNItems = i
                break
    if nItems < 0:
        raise NameError("The number of items in %s could not be extracted" % pathToFile)
    return (nItems, lineNItems)

def pathToPolyMesh(caseName, timePath=None, regionName=None):
    """
    Build the path to the polyMesh directory

    Args:
        caseName: Path to the OpenFOAM case
        timePath: Name of time folder, if the points of interest are in there
        regionName: Name of the mesh region, if any.

    Returns:
        path to polyMesh
    """
    if timePath is None and regionName is None:
        path = '%s/constant/polyMesh' % caseName
    elif timePath is None:
        path = '%s/constant/%s/polyMesh' % (caseName, regionName)

    return path


if __name__ == "__main__":
    readPoints('.')
