import operator
import sys


def getExtremumCoord(filePath, operator_str, axis):
    """
    Get the maximum coordinate of the vertices in the STL surface for a
    given axis.

    Args:
        operator_str: "min" or "max'.
        axis: coordinate axis that has to be read.
    """
    if (operator_str == 'min'):
        op = operator.lt
        extremum = 1e30
    elif (operator_str == 'max'):
        op = operator.gt
        extremum = -1e30
    else:
        raise ValueError('Error, operator must be "min" or "max".')

    with open(filePath, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            fields = [x for x in line.split()]
            if fields[0] == "vertex":
                coord = float(fields[axis+1])
                if (op(coord, extremum)):
                    extremum = coord
    return extremum


def getCenterlineExtremum(filePath, operator_str, axis, y_min, y_max):
    """
    Get the maximum coordinate of the vertices in the STL surface for a
    # given axis that have the y-coordinate between y_min and y_max.

    Args:
        operator_str: "min" or "max"
        axis: coordinate axis that has to be read.
    """
    if (operator_str == 'min'):
        op = operator.lt
        extremum = 1e30
    elif (operator_str == 'max'):
        op = operator.gt
        extremum = -1e30
    else:
        print 'Error, operator must be "min" or "max".'
        sys.exit()

    with open(filePath, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            fields = [x for x in line.split()]
            if fields[0] == "vertex":
                y = float(fields[2])
                if y > y_min and y < y_max:
                    coord = float(fields[axis+1])
                    if (op(coord, extremum)):
                        extremum = coord

    return extremum


def getLength(filePath, axis):
    max_coord = getExtremumCoord(filePath, 'max', axis)
    min_coord = getExtremumCoord(filePath, 'min', axis)
    return max_coord - min_coord


def getCenterlineLength(filePath):
    max_coord = getCenterlineExtremum(filePath, 'max', 0, -1e-7, 1e-7)
    min_coord = getCenterlineExtremum(filePath, 'min', 0, -1e-7, 1e-7)
    return max_coord - min_coord


def transformSurface(filePath, scale_factor, vector):
    """
    Transform a STL surface by scaling it and translating it.
    Write the result to stdout

    Args:
        filePath: path to STL file
        scale_factor: scaling factor
        vector: translation vector
    """
    with open(filePath, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            fields = [x for x in line.split()]
            if fields[0] == "vertex":
                x = scale_factor*float(fields[1]) + vector[0]
                y = scale_factor*float(fields[2]) + vector[1]
                z = scale_factor*float(fields[3]) + vector[2]
                print "    vertex", x, y, z
            else:
                print str(line)
