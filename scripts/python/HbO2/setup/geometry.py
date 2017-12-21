"""
Classes for geometries of computational domains for HbO2 computations.
"""

import cPickle as pickle
import json
import numpy as np
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
import warnings

from parse import readfile


class AxisymmetricGeometry(dict):
    """
    Defines the spatial dimensions of an axisymmetric geometry.
    """
    cortex_dict = {'domainLength': 300e-6,
                   'radiusTissueLeft': 26e-6,
                   'radiusTissueRight': 20e-6,
                   'radiusPlasma': 2.0e-6,
                   'radiusWall': 2.6e-6}
    glomerulus_dict = {'domainLength': 100e-6,
                      'radiusTissueLeft': 19e-6,
                      'radiusTissueRight': 13e-6,
                      'radiusPlasma': 2.0e-6,
                      'radiusWall': 2.6e-6}

    def __init__(self, geometry_dict=None):
        super(AxisymmetricGeometry, self).__init__()
        self['domainLength']      = None
        self['radiusTissueLeft']  = None
        self['radiusTissueRight'] = None
        self['radiusPlasma']      = None
        self['radiusWall']        = None
        if geometry_dict:
            for key in self:
                try:
                    self[key] = geometry_dict[key]
                except KeyError:
                    pass

    @classmethod
    def fromJson(cls, pathToJson):
        with open(pathToJson) as fp:
            jsonData = json.load(fp)
        return cls(jsonData)

    @classmethod
    def fromKeyValueFile(cls, pathToFile):
        with open(pathToFile) as fp:
            pattern = readfile.key_space_value_semicolon_pattern
            geometry_dict = readfile.readKeysFromFile(fp, pattern)
        return cls(geometry_dict)

    @classmethod
    def fromGeometryName(cls, geometryName):
        if geometryName == 'cortex':
            return cls(cls.cortex_dict)
        elif geometryName == 'glomerulus':
            return cls(cls.glomerulus_dict)

    def tissueRadius(self, x):
        L = self['domainLength']
        Ra = self['radiusTissueLeft']
        Rv = self['radiusTissueRight']
        return (1-x/L)*Ra + x/L*Rv

    def isCortex(self):
        return self.isGeometry('cortex')

    def isGlomerulus(self):
        return self.isGeometry('glomerulus')

    def isGeometry(self, geometryName):
        if hasattr(self, geometryName + '_dict'):
            geometryDict = getattr(self, geometryName + '_dict')
            if all([np.isclose(self[key], geometryDict[key]) for key in self]):
                return True
        return False


class ParallelCapillaries(object):
    """
    Define the geometry of parallel capillaries.

    All the capillaries must have the same length.
    """

    def __init__(self, dictList):
        self.geometries = []
        for dict in dictList:
            self.geometries.append(AxisymmetricGeometry(dict))
        self._checkInput()

    @classmethod
    def fromJson(cls, pathToJson):
        """
        Reads a json file with an array of objects that represent the capillary geometries
        """
        with open(pathToJson) as fp:
            jsonData = json.load(fp)
        return cls(jsonData)

    @classmethod
    def fromKeyValueFile(cls, pathToFile, nCapillaries=2):
        """
        Reads a file with key-value pairs that specify an axisymmetric geometry.

        The geometry defined in the file is duplicated n times, where n = nCapillaries.

        Args:
            pathToFile (str): path to file with key-value pairs
            nCapillaries (int): number of capillaries
        """
        geometry = AxisymmetricGeometry.fromKeyValueFile(pathToFile)
        geometries = [geometry]*nCapillaries
        return cls(geometries)

    def nCapillaries(self):
        return len(self.geometries)

    def domainLength(self):
        return self.geometries[0]['domainLength']

    def totalSliceArea(self, x):
        """Return the total area of the slice supplied by the capillaries at x.

        The capillary area is included.
        """
        return sum([np.pi*geom.tissueRadius(x)**2 for geom in self.geometries])

    def tissueSliceArea(self, x):
        """Return the area of the tissue slice supplied by the capillaries at x.

        The capillary area is not included.
        """
        return self.totalSliceArea(x) \
             - sum([np.pi*geom['radiusWall']**2 for geom in self.geometries])

    def _checkInput(self):
        if len(self.geometries) == 0:
            raise ValueError('No geometry was loaded')
        domainLengths = [geom['domainLength'] for geom in self.geometries]
        if not all(L == domainLengths[0] for L in domainLengths):
            raise ValueError('The given domain lengths are not identical')


class BoundingBoxError(Exception):
    pass

class BoundingBoxErrorAllPointsWithin(BoundingBoxError):
    pass

class BoundingBoxErrorAllPointsOutside(BoundingBoxError):
    pass

class BoundingBoxErrorNotImplemented(BoundingBoxError):
    pass


class GraphCapillaries(dict):
    """
    Representation of the graph structure of a capillary network.
    """

    required_keys = ['vertexPositions', 'segmentDiameters',
                     'adjacencyList', 'length']

    def __init__(self, graph_dict, convert_to_meters=False):
        super(GraphCapillaries, self).__init__()
        for key in graph_dict:
            self[key] = graph_dict[key]
        self._check_keys()
        if convert_to_meters:
            self._convert_to_meters()

    @classmethod
    def from_pickle(cls, path_to_pickle):
        graph_dict = pickle.load(open(path_to_pickle, 'rb'))
        # if present, rename the key for the 'tortuous' points of an edge
        if 'points' in graph_dict:
            graph_dict['edgePoints'] = graph_dict.pop('points')

    @classmethod
    def from_openfoam_dict(cls, path_to_of_dict):
        graph_dict = {}
        parsed_file = ParsedParameterFile(path_to_of_dict)
        for key, value in parsed_file.getValueDict().iteritems():
            graph_dict[key] = value
        graph_dict['vertexPositions'] = [elem.vals for elem in graph_dict['vertexPositions']]
        if 'edgePoints' in graph_dict:
            graph_dict['edgePoints'] = [np.array([elem.vals for elem in point_list])
                                                 for point_list in graph_dict['edgePoints']]
        # compute length
        graph_dict['length'] = graph_edge_lengths(graph_dict)
        return cls(graph_dict, convert_to_meters=False)

    def n_edges(self):
        return len(self['adjacencyList'])

    def edge_length(self, ei):
        """
        Return the length of the edge with index ei

        Args:
            ei (int): edge index

        Returns:
            float, edge length
        """
        return self['length'][ei]

    def edge_diameter(self, ei):
        return np.mean(self['segmentDiameters'][ei])

    def edge_radius(self, ei):
        return 0.5*self.edge_diameter(ei)

    def edge_coordinate_interval_in_box(self, ei, box_min, box_max):
        """
        Return the bounds of the curvilinear coordinates on an edge with index ei
        that are contained within a bounding box given by min_point and max_point.

        Args:
            ei (int): edge index
            box_min (np.ndarray): minimum point of the bounding box
            box_max (np.ndarray): maximum point of the bounding box

        Returns:
            2-tuple of float, with ordered coordinates.

        Raises:
            BoundingBoxError if the edge is outside the bounding box or if it crosses it multiple times.
        """
        try:
            edge_points = self['edgePoints'][ei]
        except KeyError:
            edge_points = np.array([self['vertexPositions'][self['adjacencyList'][ei][0]],
                                    self['vertexPositions'][self['adjacencyList'][ei][1]]])
        n_points = edge_points.shape[0]
        in_box_mask = [is_in_box(edge_points[i,:], box_min, box_max) for i in range(n_points)]
        if all(in_box_mask):
            return (0, self['length'][ei])
        elif not any(in_box_mask):
            raise BoundingBoxErrorAllPointsOutside('All segment points are outside the bounding box')
        else:
            cutting_segment_ids = [i for i in range(n_points - 1) if in_box_mask[i] != in_box_mask[i+1]]
            if len(cutting_segment_ids) >= 2:
                warnings.warn("The edge {:d} cuts the bounding box more than once.\n".format(ei),
                              UserWarning)
            if in_box_mask[0]:
                length = 0
                for i in range(n_points - 1):
                    if in_box_mask[i] and in_box_mask[i+1]:
                        length += np.linalg.norm(edge_points[i+1] - edge_points[i])
                    elif in_box_mask[i] and not in_box_mask[i+1]:
                        p = segment_intersection_with_box(edge_points[i], edge_points[i+1],
                                                          box_min, box_max)
                        length += np.linalg.norm(p - edge_points[i])
                return (0.0, length)
            elif in_box_mask[-1]:
                length = 0
                entering_coord = -1
                for i in range(n_points - 1):
                    if not in_box_mask[i] and not in_box_mask[i+1]:
                        length += np.linalg.norm(edge_points[i+1] - edge_points[i])
                    if not in_box_mask[i] and in_box_mask[i+1]:
                        p = segment_intersection_with_box(edge_points[i], edge_points[i+1],
                                                          box_min, box_max)
                        entering_coord = length + np.linalg.norm(p - edge_points[i])
                return entering_coord, self['length'][ei]

            raise BoundingBoxErrorNotImplemented("This case is not implemented yet")

    def edge_coordinate_intervals_in_box(self, ei, box_min, box_max, buffer_width=0.0):
        """
        Return the bounds of the curvilinear coordinates on an edge with index ei
        that are contained within a bounding box given by min_point and max_point.

        A buffer zone specified by its width can be used so that an edge that leaves the box
        but stays within the buffer zone (and possibly reenters the box) is not considered
        to leave the box. This is useful if the considered edges are actually tubes with a
        certain radius.

        Args:
            ei (int): edge index
            box_min (np.ndarray): minimum point of the bounding box
            box_max (np.ndarray): maximum point of the bounding box
            buffer_width (float): width of the buffer zone

        Returns:
            out (list): list of 2-tuples of floats with ordered coordinate intervals

        Raises:
            BoundingBoxError if the edge is outside the bounding box
        """
        try:
            edge_points = self['edgePoints'][ei]
        except KeyError:
            edge_points = np.array([self['vertexPositions'][self['adjacencyList'][ei][0]],
                                    self['vertexPositions'][self['adjacencyList'][ei][1]]])
        n_points = edge_points.shape[0]
        in_box_mask = [is_in_box(edge_points[i, :], box_min, box_max) for i in range(n_points)]
        buffer_box_min = box_min - buffer_width
        buffer_box_max = box_max + buffer_width
        in_buffer_mask = [is_in_box(edge_points[i, :], buffer_box_min, buffer_box_max)
                          and not in_box_mask[i] for i in range(n_points)]
        if all(in_box_mask):
            return [(0, self['length'][ei])]
        else:
            intervals = []
            cum_length = 0.0
            s_min = 0.0
            is_inside = in_box_mask[0]   # whether the i-th point is inside the bounding box
            for segment_i in range(n_points - 1):
                is_staying = is_staying_segment(segment_i, in_box_mask, in_buffer_mask)
                if is_staying:
                    is_inside = True
                else:
                    try:
                        intersection_points = segment_intersections_with_box(edge_points[segment_i],
                                                                             edge_points[segment_i+1],
                                                                             box_min, box_max)
                        n_inters = 1 if intersection_points.shape == (3,) else 2
                        is_entering = is_entering_segment(segment_i, in_box_mask, in_buffer_mask)
                        is_leaving = is_leaving_segment(segment_i, in_box_mask, in_buffer_mask)
                        if is_leaving:
                            s_max = cum_length + np.linalg.norm(intersection_points - edge_points[segment_i])
                            is_inside = False
                            intervals.append((s_min, s_max))
                        elif is_entering:
                            s_min = cum_length + np.linalg.norm(intersection_points - edge_points[segment_i])
                            is_inside = True
                        elif n_inters == 2:
                            s_min = cum_length + np.linalg.norm(intersection_points[0, :] - edge_points[segment_i])
                            s_max = cum_length + np.linalg.norm(intersection_points[1, :] - edge_points[segment_i])
                            is_inside = False
                            intervals.append((s_min, s_max))
                    except BoundingBoxErrorAllPointsOutside:  # whole segment is outside
                        is_inside = False
                cum_length += np.linalg.norm(edge_points[segment_i+1] - edge_points[segment_i])
            if is_inside:
                s_max = cum_length
                intervals.append((s_min, s_max))
        return intervals

    def _check_keys(self):
        if not set(self.keys()).issuperset(set(self.required_keys)):
            print "Required keys for GraphCapillaries instance:"
            print self.required_keys
            print "Present keys:"
            print self.keys()
            raise RuntimeError('A required key is missing.')

    def _convert_to_meters(self):
        """
        Convert the dimensions in the graph to meters if they are given in micrometers.
        """
        self['length'] = [1e-6*l for l in self['length']]
        self['vertexPositions'] = [1e-6*p for p in self['vertexPositions']]
        self['segmentDiameters'] = [1e-6*np.asarray(d) for d in self['segmentDiameters']]
        if 'edgePoints' in self:
            self['edgePoints'] = [1e-6*edge_pts for edge_pts in self['edgePoints']]


def is_in_box(p, box_min, box_max):
    """
    Return whether the given point is in the bounding box given by box_min and box_max.

    Args:
        p (np.ndarray): point to test
        box_min (np.ndarray): minimum point of the bounding box
        box_max (np.ndarray): maximum point of the bounding box

    Returns:
        bool
    """
    return all(p >= box_min) and all(p <= box_max)


def segment_intersection_with_box(p1, p2, box_min, box_max):
    """
    Return the intersection point of the segment between p1 and p2 with the bounding box
    given by box_min and box_max. One point should be inside the bounding box and the
    other outside it.

    Uses a binary search, so the solution is not exact.

    Args:
        p1 (np.ndarray): segment end point
        p2 (np.ndarray): segment end point
        box_min (np.ndarray): minimum point of the bounding box
        box_max (np.ndarray): maximum point of the bounding box

    Returns:
        np.ndarray, intersection point
    """
    if not is_in_box(p1, box_min, box_max) and not is_in_box(p2, box_min, box_max):
        raise BoundingBoxErrorAllPointsOutside('Both segment end points are outside the bounding box')
    elif is_in_box(p1, box_min, box_max) and is_in_box(p2, box_min, box_max):
        raise BoundingBoxErrorAllPointsWithin('Both segment end points are within the bounding box')

    q_in  = p1 if is_in_box(p1, box_min, box_max) else p2
    q_out = p2 if is_in_box(p1, box_min, box_max) else p1
    dist = np.inf
    tol = 1e-8 * np.max(np.abs(box_max - box_min))
    while dist > tol:
        mid_point = 0.5*(q_in + q_out)
        if is_in_box(mid_point, box_min, box_max):
            q_in = mid_point
        else:
            q_out = mid_point
        dist = np.linalg.norm(q_in - q_out)
    return 0.5*(q_in + q_out)


def is_leaving_segment(i, in_box_mask, in_buffer_mask):
    n_points = len(in_box_mask)
    if in_box_mask[i] and not in_box_mask[i+1]:
        if not in_buffer_mask[i+1]:
            return True
        else:
            try:
                id_leaving_buffer = [in_buffer_mask[j] and not in_buffer_mask[j+1]
                                     for j in range(i + 1, n_points - 1)].index(True) + i + 1
                if in_box_mask[id_leaving_buffer+1]:
                    return False
                else:
                    return True
            except ValueError:
                return False
    else:
        return False


def is_entering_segment(i, in_box_mask, in_buffer_mask):
    if not in_box_mask[i] and in_box_mask[i+1]:
        if not in_buffer_mask[i]:
            return True
        else:
            try:
                id_entering_buffer = [not in_buffer_mask[j] and in_buffer_mask[j+1]
                                     for j in reversed(range(i))].index(True)*(-1) + i - 1
                if not in_box_mask[id_entering_buffer]:
                    return True
                else:
                    return False
            except ValueError:
                return True
    else:
        return False


def is_staying_segment(i, in_box_mask, in_buffer_mask):
    is_outside_mask = [not in_box and not in_buffer
                       for in_box, in_buffer in zip(in_box_mask, in_buffer_mask)]
    if is_outside_mask[i] or is_outside_mask[i+1]:
        return False
    if in_box_mask[i] and in_box_mask[i+1]:
        return True
    elif in_box_mask[i] and in_buffer_mask[i+1] \
        and not is_leaving_segment(i, in_box_mask, in_buffer_mask):
        return True
    elif in_buffer_mask[i] and in_box_mask[i+1] \
            and not is_entering_segment(i, in_box_mask, in_buffer_mask):
        return True
    else:
        try:
            id_entering_buffer = [not in_buffer_mask[j] and in_buffer_mask[j+1]
                                  for j in reversed(range(i))].index(True)*(-1) + i - 1
            inside_before = in_box_mask[id_entering_buffer]
        except ValueError:
            inside_before = True
        try:
            n_points = len(in_box_mask)
            id_leaving_buffer = [in_buffer_mask[j] and not in_buffer_mask[j+1]
                                 for j in range(i + 1, n_points - 1)].index(True) + i + 1
            inside_after = in_box_mask[id_leaving_buffer]
        except ValueError:
            inside_after = True
        return inside_before and inside_after


def segment_intersections_with_box(p1, p2, box_min, box_max, min_length=1e-6):
    """
    Return the intersection points of the segment between p1 and p2 with the bounding box
    given by box_min and box_max. One point should be inside the bounding box and the
    other outside it.

    Uses a binary search, so the solution is not exact.

    Args:
        p1 (np.ndarray): segment end point
        p2 (np.ndarray): segment end point
        box_min (np.ndarray): minimum point of the bounding box
        box_max (np.ndarray): maximum point of the bounding box
        min_length (float): minimum segment length resolved by the algorithm

    Returns:
        np.ndarray, intersection points with shape (n, 3), where n is the number of
        intersection points
    """
    if is_in_box(p1, box_min, box_max) and is_in_box(p2, box_min, box_max):
        raise BoundingBoxErrorAllPointsWithin('Both segment end points are within the bounding box')
    # if the end points do not cross the planes that define the bounding box, raise an error
    if any([   (p1[i] < box_min[i] and p2[i] < box_min[i])
            or (p1[i] > box_max[i] and p2[i] > box_max[i]) for i in range(3)]):
        raise BoundingBoxErrorAllPointsOutside

    try:
        return segment_intersection_with_box(p1, p2, box_min, box_max)
    except BoundingBoxErrorAllPointsOutside:
        # termination condition for the recursion
        if np.linalg.norm(p1 - p2) < min_length:
            raise BoundingBoxErrorAllPointsOutside
        # search for an intersection point on both sides of the midpoint
        mid_point = 0.5*(p1+p2)
        q1 = np.array(())
        q2 = np.array(())
        q1_found = False
        q2_found = False
        try:
            q1 = segment_intersections_with_box(p1, mid_point, box_min, box_max, min_length)
            q1_found = True
        except BoundingBoxError:
            pass
        try:
            q2 = segment_intersections_with_box(mid_point, p2, box_min, box_max, min_length)
            q2_found = True
        except BoundingBoxError:
            pass

        if q1_found and q2_found:
            return np.stack([q1, q2], axis=0)
        elif q1_found and not q2_found:
            return q1
        elif not q1_found and q2_found:
            return q2
        else:
            raise BoundingBoxErrorAllPointsOutside


def graph_edge_lengths(graph_dict):
    if 'edgePoints' in graph_dict:
        return [sum([np.linalg.norm(np.array(point_list[i+1]) - np.array(point_list[i]))
                     for i in range(len(point_list) - 1)])
                for point_list in graph_dict['edgePoints']]
    else:
        positions = np.asarray(graph_dict['vertexPositions'])
        return [np.linalg.norm(positions[e[1],:] - positions[e[0],:])
                    for e in graph_dict['adjacencyList']]