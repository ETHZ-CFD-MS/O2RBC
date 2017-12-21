#!/usr/bin/env python
"""
Unit test for HbO2/setup/geometry.py
"""

import numpy as np
import unittest

from HbO2.setup import geometry


class IsInBox(unittest.TestCase):
    def testIsInBox(self):
        p1 = np.array([0.5, 0.5, 0.5])
        p2 = np.array([1.5, 0.5, 0.5])
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        self.assertTrue(geometry.is_in_box(p1, box_min, box_max))
        self.assertFalse(geometry.is_in_box(p2, box_min, box_max))


class SegmentIntersectionWithBox(unittest.TestCase):
    def testIntersectionPoint(self):
        p1 = np.array([-0.25, 0.25, 0.25])
        p2 = np.array([0.25,  0.75, 0.75])
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        self.assertAlmostEqual(
            np.linalg.norm(
                geometry.segment_intersection_with_box(p1, p2, box_min, box_max) -
                np.array([0.0, 0.5, 0.5])),
            0.0,
            places=8)
        q1 = np.array([0.25,   0.25,  0.25])
        q2 = np.array([-0.25, -0.25, -0.25])
        self.assertAlmostEqual(
            np.linalg.norm(
                geometry.segment_intersection_with_box(q1, q2, box_min, box_max) -
                np.array([0.0, 0.0, 0.0])),
            0.0,
            places=8)

    def testBothPointsOutside(self):
        p1 = np.array([-0.5, 0.5, 0.5])
        p2 = np.array([-0.25, 0.5, 0.5])
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        self.assertRaises(geometry.BoundingBoxErrorAllPointsOutside,
                          geometry.segment_intersection_with_box,
                          p1, p2, box_min, box_max)

    def testBothPointsInside(self):
        p1 = np.array([0.5, 0.5, 0.5])
        p2 = np.array([0.25, 0.5, 0.5])
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        self.assertRaises(geometry.BoundingBoxErrorAllPointsWithin,
                          geometry.segment_intersection_with_box,
                          p1, p2, box_min, box_max)


class SegmentIntersectionsWithBox(unittest.TestCase):

    def testOneIntersectionPoints(self):
        p1 = np.array([-0.25, 0.5, 0.5])
        p2 = np.array([0.25,  0.5, 0.5])
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        intersection_points = geometry.segment_intersections_with_box(p1, p2, box_min, box_max)
        self.assertAlmostEqual(
            np.linalg.norm(intersection_points - np.array([0.0, 0.5, 0.5])),
            0.0,
            places=8)

    def testTwoIntersectionPoints(self):
        p1 = np.array([-0.25, 0.5, 0.5])
        p2 = np.array([1.25,  0.5, 0.5])
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        intersection_points = geometry.segment_intersections_with_box(p1, p2, box_min, box_max)
        self.assertAlmostEqual(
            np.linalg.norm(intersection_points[0, :] - np.array([0.0, 0.5, 0.5])),
            0.0,
            places=8)
        self.assertAlmostEqual(
            np.linalg.norm(intersection_points[1, :] - np.array([1.0, 0.5, 0.5])),
            0.0,
            places=8)
        p1 = np.array([-1e-4 + 1e-6, 1e-4, 0.5])
        p2 = np.array([ 1e-3, -1e-3 + 1e-6, 0.5])
        intersection_points = geometry.segment_intersections_with_box(p1, p2, box_min, box_max)
        self.assertAlmostEqual(
            np.linalg.norm(intersection_points[0, :] - np.array([0.0, 1e-6, 0.5])),
            0.0,
            places=8)
        self.assertAlmostEqual(
            np.linalg.norm(intersection_points[1, :] - np.array([1e-6, 0.0, 0.5])),
            0.0,
            places=8)

    def testAllPointsOutside(self):
        # simple case where the points do not cross any of the planes that define the bounding box
        p1 = np.array([-0.5, 0.5, 0.5])
        p2 = np.array([-0.4, 0.5, 0.5])
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        self.assertRaises(geometry.BoundingBoxErrorAllPointsOutside,
                          geometry.segment_intersections_with_box,
                          p1, p2, box_min, box_max)
        p1 = np.array([-0.2, 0.1, -0.1])
        p2 = np.array([ 0.1, -0.2, 0.1])
        self.assertRaises(geometry.BoundingBoxErrorAllPointsOutside,
                          geometry.segment_intersections_with_box,
                          p1, p2, box_min, box_max)

class EdgeCoordinateIntervalInBox(unittest.TestCase):

    def testOneEnteringEdge(self):
        graph_dict = {'adjacencyList': [(0, 1)],
                      'segmentDiameters': [1],
                      'vertexPositions': [np.array([0.5, -0.25, 0.5]), np.array([0.5, 0.5, 0.5])],
                      'length': [0.75]}
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        graph_geometry = geometry.GraphCapillaries(graph_dict)
        coords = graph_geometry.edge_coordinate_interval_in_box(0, box_min, box_max)
        self.assertAlmostEqual(coords[0], 0.25)
        self.assertAlmostEqual(coords[1], 0.75)

    def testOneLeavingEdge(self):
        graph_dict = {'adjacencyList': [(0, 1)],
                      'segmentDiameters': [1],
                      'vertexPositions': [np.array([0.5, 0.25, 0.5]), np.array([0.5, -0.25, 0.5])],
                      'length': [0.5]}
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        graph_geometry = geometry.GraphCapillaries(graph_dict)
        coords = graph_geometry.edge_coordinate_interval_in_box(0, box_min, box_max)
        self.assertAlmostEqual(coords[0], 0)
        self.assertAlmostEqual(coords[1], 0.25)

    def testMultiSegmentEnteringEdge(self):
        graph_dict = {'adjacencyList': [(0, 1)],
                      'segmentDiameters': [1],
                      'vertexPositions': [np.array([0.5, -0.25, 0.5]), np.array([0.5, 0.5, 0.5])],
                      'points': [np.array([[0.5, -0.25, 0.5],[0.5, -0.15, 0.5],
                                           [0.5, 0.25, 0.5], [0.5, 0.5, 0.5]])],
                      'length': [0.75]}
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        graph_geometry = geometry.GraphCapillaries(graph_dict)
        coords = graph_geometry.edge_coordinate_interval_in_box(0, box_min, box_max)
        self.assertAlmostEqual(coords[0], 0.25)
        self.assertAlmostEqual(coords[1], 0.75)

    def testMultiSegmentLeavingEdge(self):
        graph_dict = {'adjacencyList': [(0, 1)],
                      'segmentDiameters': [1],
                      'vertexPositions': [np.array([0.5, 0.5, 0.5]), np.array([0.5, -0.25, 0.5])],
                      'points': [np.array([[0.5, 0.5, 0.5], [0.5, 0.25, 0.5], [0.5, -0.25, 0.5]])],
                      'length': [0.75]}
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        graph_geometry = geometry.GraphCapillaries(graph_dict)
        coords = graph_geometry.edge_coordinate_interval_in_box(0, box_min, box_max)
        self.assertAlmostEqual(coords[0], 0)
        self.assertAlmostEqual(coords[1], 0.5)

    def testContainedEdge(self):
        graph_dict = {'adjacencyList': [(0, 1)],
                      'segmentDiameters': [1],
                      'vertexPositions': [np.array([0.5, 0.25, 0.5]), np.array([0.5, 0.75, 0.5])],
                      'length': [0.5]}
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        graph_geometry = geometry.GraphCapillaries(graph_dict)
        coords = graph_geometry.edge_coordinate_interval_in_box(0, box_min, box_max)
        self.assertAlmostEqual(coords[0], 0)
        self.assertAlmostEqual(coords[1], 0.5)

    def testExteriorEdge(self):
        graph_dict = {'adjacencyList': [(0, 1)],
                      'segmentDiameters': [1],
                      'vertexPositions': [np.array([-1.5, 0.5, 0.5]), np.array([-0.5, 0.5, 0.5])],
                      'length': [1]}
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        graph_geometry = geometry.GraphCapillaries(graph_dict)
        self.assertRaises(geometry.BoundingBoxError, graph_geometry.edge_coordinate_interval_in_box,
                          0, box_min, box_max)

class EdgeCoordinateIntervalsInBox(unittest.TestCase):

    def testOneEnteringEdge(self):
        graph_dict = {'adjacencyList': [(0, 1)],
                      'segmentDiameters': [1],
                      'vertexPositions': [np.array([0.5, -0.25, 0.5]), np.array([0.5, 0.5, 0.5])],
                      'length': [0.75]}
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        graph_geometry = geometry.GraphCapillaries(graph_dict)
        intervals = graph_geometry.edge_coordinate_intervals_in_box(0, box_min, box_max)
        self.assertAlmostEqual(intervals[0][0], 0.25)
        self.assertAlmostEqual(intervals[0][1], 0.75)

    def testOneLeavingEdge(self):
        graph_dict = {'adjacencyList': [(0, 1)],
                      'segmentDiameters': [1],
                      'vertexPositions': [np.array([0.5, 0.25, 0.5]), np.array([0.5, -0.25, 0.5])],
                      'length': [0.5]}
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        graph_geometry = geometry.GraphCapillaries(graph_dict)
        intervals = graph_geometry.edge_coordinate_intervals_in_box(0, box_min, box_max)
        self.assertAlmostEqual(intervals[0][0], 0)
        self.assertAlmostEqual(intervals[0][1], 0.25)

    def testEnteringLeavingEdge(self):
        graph_dict = {'adjacencyList': [(0, 1)],
                      'segmentDiameters': [1],
                      'vertexPositions': [np.array([-0.1, 0.5, 0.5]), np.array([1.1, 0.5, 0.5])],
                      'length': [1.2]}
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        graph_geometry = geometry.GraphCapillaries(graph_dict)
        intervals = graph_geometry.edge_coordinate_intervals_in_box(0, box_min, box_max)
        self.assertAlmostEqual(intervals[0][0], 0.1)
        self.assertAlmostEqual(intervals[0][1], 1.1)

    def testLeavingEnteringEdge(self):
        graph_dict = {'adjacencyList': [(0, 1)],
                      'segmentDiameters': [1],
                      'vertexPositions': [np.array([0.25, 0.5, 0.5]), np.array([0.75, 0.5, 0.5])],
                      'points': [np.array([[0.25, 0.5, 0.5], [0.25, 1.5, 0.5],
                                           [0.75, 1.5, 0.5], [0.75, 0.5, 0.5]])],
                      'length': [2.5]}
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        graph_geometry = geometry.GraphCapillaries(graph_dict)
        intervals = graph_geometry.edge_coordinate_intervals_in_box(0, box_min, box_max)
        self.assertAlmostEqual(intervals[0][0], 0.0)
        self.assertAlmostEqual(intervals[0][1], 0.5)
        self.assertAlmostEqual(intervals[1][0], 2.0)
        self.assertAlmostEqual(intervals[1][1], 2.5)

    def testMultipleEnteringLeavingEdge(self):
        graph_dict = {'adjacencyList': [(0, 1)],
                      'segmentDiameters': [1],
                      'vertexPositions': [np.array([-0.25, 0.5, 0.5]), np.array([0.75, 0.5, 0.5])],
                      'points': [np.array([[-0.25, 0.5, 0.5], [-0.1, 0.5, 0.5],
                                           [0.25, 0.5, 0.5], [0.25, 1.5, 0.5],
                                           [0.5, 1.5, 0.5], [0.75, 1.5, 0.5],
                                           [0.75, 0.5, 0.5]])],
                      'length': [3.0]}
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        graph_geometry = geometry.GraphCapillaries(graph_dict)
        intervals = graph_geometry.edge_coordinate_intervals_in_box(0, box_min, box_max)
        self.assertAlmostEqual(intervals[0][0], 0.25)
        self.assertAlmostEqual(intervals[0][1], 1.0)
        self.assertAlmostEqual(intervals[1][0], 2.5)
        self.assertAlmostEqual(intervals[1][1], 3.0)

    def testStayingEdgeWidthBuffer(self):
        graph_dict = {'adjacencyList': [(0, 1)],
                      'segmentDiameters': [1],
                      'vertexPositions': [np.array([0.25, 0.5, 0.5]), np.array([0.75, 0.5, 0.5])],
                      'points': [np.array([[0.25, 0.5, 0.5], [0.25, 1.1, 0.5],
                                           [0.75, 1.1, 0.5], [0.75, 0.5, 0.5]])],
                      'length': [1.7]}
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([1.0, 1.0, 1.0])
        graph_geometry = geometry.GraphCapillaries(graph_dict)
        intervals = graph_geometry.edge_coordinate_intervals_in_box(0, box_min, box_max,
                                                                    buffer_width=0.2)
        self.assertAlmostEqual(intervals[0][0], 0)
        self.assertAlmostEqual(intervals[0][1], 1.7)

        intervals2 = graph_geometry.edge_coordinate_intervals_in_box(0, box_min, box_max,
                                                                     buffer_width=0.05)
        self.assertAlmostEqual(intervals2[0][0], 0)
        self.assertAlmostEqual(intervals2[0][1], 0.5)
        self.assertAlmostEqual(intervals2[1][0], 1.2)
        self.assertAlmostEqual(intervals2[1][1], 1.7)


if __name__ == '__main__':
    unittest.main()