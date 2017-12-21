import functools
import json
import numpy as np
import warnings

from utilities.decorators import lazy_function, lazy_property, remove_lazy_properties


class RBCDataPostProcessor(object):
    "Provides methods for the postprocessing of RBC data"

    def __init__(self, rbc_data):
        """
        Args:
            rbc_data (data.rbcdata.RBCData): RBC data
        """
        self.rbc_data = rbc_data

    def fieldAverage(self, fieldName, x, nAverage=10):
        """Return the averaged field with name field_name at position x,
        averaged over the last nAverage paths"""
        field = self.fieldOnLastPaths(fieldName, x, nAverage)
        averaging_axis = min(1, field.ndim-1)
        return np.mean(field, axis=averaging_axis)

    def fieldAverageAlternating(self, fieldName, x, nGroups, nAverage=10):
        """Return field at x averaged over groups of alternating values.
        nAverage is the number of RBC paths per group.
        The returned array has shape (len(x), nGroups)."""
        nPath = nGroups*nAverage
        values = self.fieldOnLastPaths(fieldName, x, nPath)
        averages = np.zeros( (len(x), nGroups) )
        for i in range(nGroups):
            indices = range(i, values.shape[1], nGroups)   # use the actual number
                                                           # of paths rather than nPath
            averages[:,i] = np.mean(values[:,indices], axis=1)
        return averages

    def fieldStd(self, fieldName, x, nAverage=20):
        """Return the standard deviation of the given field at position x,
        computed over the nAverage last paths."""
        field = self.fieldOnLastPaths(fieldName, x, nAverage)
        averaging_axis = min(1, field.ndim-1)
        return np.std(field, axis=averaging_axis)

    def fieldOnLastPaths(self, fieldName, x, nPath):
        """Returns a numpy array with the field values at x for the last
        nPath RBC paths. The array shape is (len(x), nPath)."""
        if nPath > self.rbc_data.nPath():
            warnings.warn(
                """The required number of paths (%i) is larger than the number
                of available paths (%i). Using all available paths.
                """ % (nPath, self.rbc_data.nPath()), UserWarning)
            nPath = self.rbc_data.nPath()
        pathIndices = range(-nPath, 0)
        f = functools.partial(self.rbc_data.fieldValue, fieldName, x)
        return np.asarray(map(f, pathIndices)).transpose()


class RBCDataGraphPostProcessor(RBCDataPostProcessor):
    """
    Postprocess RBC data for O2 simulations in graphs.
    """

    def __init__(self, rbc_data):
        super(RBCDataGraphPostProcessor, self).__init__(rbc_data)
        self.rbc_path_analyzer = RBCPathAnalyzer(rbc_data)

    def fieldOnEdge(self, fieldName, s, eI, nPath):
        """
        Return the n last field values on a given edge for given coordinates

        Args:
            fieldName (str): field name
            s (numpy.ndarray): curvilinear coordinate on edge
            eI (int): edge index
            nPath (int): number of RBCs for which to return values

        Returns:
            numpy array (len(s) x n) with field values
        """
        RBCIndices = self.rbc_path_analyzer.last_complete_path_indices_on_edge(eI, nPath)
        f = functools.partial(self.rbc_data.fieldValueOnEdge, fieldName, s, eI)
        return np.asarray(map(f, RBCIndices)).transpose()

    def fieldAverageOnEdge(self, fieldName, s, eI, nAverage=10):
        """
        Compute a field average at the given positions on a given edge.

        Args:
            fieldName (str): field name
            s (numpy.ndarray): curvilinear coordinate on edge
            eI (int): edge index
            nAverage (int): number of RBCs used for averaging

        Returns:
            np.ndarray, averaged field value at the given positions,
        """
        field = self.fieldOnEdge(fieldName, s, eI, nAverage)
        averaging_axis = min(1, field.ndim-1)
        return np.nanmean(field, axis=averaging_axis)

    def fieldStdOnEdge(self, fieldName, s, eI, nAverage=10):
        """
        Compute the field standard deviation at the given positions on a given edge.

        Args:
            fieldName (str): field name
            s (numpy.ndarray): curvilinear coordinate on edge
            eI (int): edge index
            nAverage (int): number of RBCs used for averaging

        Returns:
            np.ndarray, field standard deviation at the given positions
        """
        field = self.fieldOnEdge(fieldName, s, eI, nAverage)
        averaging_axis = min(1, field.ndim-1)
        return np.nanstd(field, axis=averaging_axis)

    def sCoordIntervalWithDefinedValues(self, fieldName, s, eI, nPath=10):
        """
        Return the s-coordinate interval for which all s-values are defined for the
        nPath last paths.

        Args:
            fieldName (str): field name
            s (numpy.ndarray): curvilinear coordinate on edge
            eI (int): edge index
            nPath (int): number of tested RBCs

        Returns:
            float length-2 tuple, interval bounds
        """
        y = np.mean(self.fieldOnEdge(fieldName, s, eI, nPath), axis=1)
        numeric_indices = np.invert(np.isnan(y))
        return (np.min(s[numeric_indices]), np.max(s[numeric_indices]))

    def timesOnEdgeCoordinate(self, s, eI, nPath):
        """
        Compute the times at which a RBC passes at the given coordinate for the last nPath RBC paths

        Args:
            s (float): curvilinear coordinate
            eI (int): edge index
            nPath (int): number of RBCs

        Returns:
            float list with time values
        """
        RBCIndices = self.rbc_path_analyzer.last_complete_path_indices_on_edge(eI, nPath)
        return [self.rbc_path_analyzer.edge_coordinate_to_time(s, eI, pathI) for pathI in RBCIndices]

    def edgeIdsWithPaths(self, nPathMin=1):
        """
        Return a list of edge indices which have a minimal number of paths crossing them.

        Args:
            nPathMin (int): minimal number of paths required so that an edge be included

        Returns:
            list of int, edge indices
        """
        return self.rbc_path_analyzer.edge_ids_with_paths(n_path_min=nPathMin)

    def timeToPreviousRBC(self, eI, nPath):
        """
        Return an array with contains the passing time difference between the last nPath RBCs
        in edge eI.

        Args:
            eI (int): edge index
            nPath (int): number of paths

        Returns:
            np.ndarray with length nPath - 1
        """
        RBCIndices = self.rbc_path_analyzer.last_complete_path_indices_on_edge(eI, nPath)
        # choose an s-coordinate that is common to all selected paths by computing the
        # largest coordinate interval contains in all paths and taking its midpoint
        s_lower_bound = max([np.min(self.rbc_data.pathCentersOnEdge(rbc_i, eI))
                             for rbc_i in RBCIndices])
        s_upper_bound = min([np.max(self.rbc_data.pathCentersOnEdge(rbc_i, eI))
                             for rbc_i in RBCIndices])
        if s_upper_bound < s_lower_bound:
            raise RuntimeError('No coordinate common to all selected paths on edge {:d}'.format(eI))
        s = 0.5*(s_lower_bound + s_upper_bound)
        times = self.timesOnEdgeCoordinate(s, eI, nPath)
        return np.diff(times)

    def removePathsWithUndefinedFlowDirection(self):
        """
        Remove the RBC paths that do not have a defined flow direction from the RBC data.
        """
        rbc_names_to_remove = set()
        for path_i in range(self.rbc_data.nPath()):
            for ei in self.rbc_data.edgesOnPath(path_i):
                if not self.rbc_path_analyzer.defined_path_flow_direction_on_edge(path_i, ei):
                    rbc_names_to_remove.add(self.rbc_data.pathIndexToRBCName(path_i))
        if len(rbc_names_to_remove) > 0:
            rbc_name_string = ', '.join([rbc_name for rbc_name in rbc_names_to_remove])
            warnings.warn('Removing the following RBCs with undefined flow direction: {:s}'
                          .format(rbc_name_string), UserWarning)
            for rbc_name in rbc_names_to_remove:
                self.rbc_data.deletePath(rbc_name)
            remove_lazy_properties(self.rbc_path_analyzer)


class FlowReversalError(RuntimeError):
    pass


class CoordinateNotFoundError(RuntimeError):
    pass


class RBCPathAnalyzer(object):
    """
    Analyze RBC paths produced by an Hb-O2 simulation in a graph.

    In contrast to a purely topological graph analysis (which could be done within vgm),
    this class analyzes the graph based on the actual RBCs that were used in the OpenFOAM
    simulation. Therefore, a connecting edge in the underlying VGM graph can be an inlet
    edge in the OpenFOAM simulation.
    """

    def __init__(self, rbc_data):
        self.rbc_data = rbc_data

    @lazy_function
    def positive_flow(self, ei, tol=1e-15):
        """
        Return True if the RBCs flow in increasing s-coordinate direction in edge ei.
        The tolerance allows for apparent flow rehearsal due to stalling RBCs.

        Args:
            ei (int): edge index
            tol (float): tolerance for testing flow reversal

        Returns:
            bool

        Raises:
            FlowReversalError if flow direction is inconsistent.
        """
        path_ids = self.rbc_data.edgeIntersectingPathIds(ei)
        positive_flow = {}
        for path_i in path_ids:
            positive_flow[path_i] = self.positive_path_flow_on_edge(path_i, ei, tol)
        if all(positive_flow.values()):
            return True
        elif not any(positive_flow.values()):
            return False
        else:
            raise FlowReversalError('Flow reversal in edge {:d}'.format(ei))

    def positive_path_flow_on_edge(self, path_i, ei, tol=1e-15):
        """
        Return True if the RBCs flow in increasing s-coordinate direction in edge ei.
        The tolerance allows for apparent flow rehearsal due to stalling RBCs.

        Args:
            path_i (int): path index
            ei (int): edge index
            tol (float): tolerance for testing flow reversal

        Returns:
            bool

        Raises:
            FlowReversalError if flow direction is inconsistent.
        """
        s = self.rbc_data.pathCentersOnEdge(path_i, ei)
        if len(s) >= 2:
            if (np.diff(s) >= -tol).all():
                return True
            elif (np.diff(s) <= tol).all():
                return False
            else:
                rbc_name = self.rbc_data.pathIndexToRBCName(path_i)
                warnings.warn('Flow reversal for {:s} in edge {:d}'.format(rbc_name, ei))
                s_low, s_high = self.cropped_s_coord_interval_with_paths(ei)
                if s[0] <= s_low and s_high <= s[-1]:
                    return True
                elif s[-1] <= s_low and s_high <= s[0]:
                    return False
                else:
                    raise FlowReversalError('Flow reversal for {:s} in edge {:d}'.
                                            format(rbc_name, ei))

    def defined_path_flow_direction_on_edge(self, path_i, ei):
        """
        Return whether the flow direction of the given path on the given edge is defined.

        Args:
            path_i (int): path index
            ei (int): edge index

        Returns:
            bool
        """
        try:
            self.positive_path_flow_on_edge(path_i, ei)
            return True
        except FlowReversalError:
            return False

    def inlet_edges(self):
        """
        Return the inlet edge indices (list of integers).
        """
        eids = self.edge_ids_with_paths()
        # remove all edges that are a successor of another edge
        for ei_succ in self.edge_successors.values():
            eids = [ei for ei in eids if ei not in ei_succ]
        return eids

    def outlet_edges(self):
        """
        Return the outlet edge indices (list of integers).
        """
        eids = self.edge_ids_with_paths()
        # remove all edges that are a predecessor of another edge
        for ei_pred in self.edge_predecessors.values():
            eids = [ei for ei in eids if ei not in ei_pred]
        return eids

    def post_converging_bifurcation_edges(self):
        """
        Return the edge indices which are after a converging bifurcation (list of integers).
        """
        return [ei for ei in self.edge_ids_with_paths() if len(self.edge_predecessors[ei]) > 1]

    def post_diverging_bifurcation_edges(self):
        """
        Return the edge indices which are after a diverging bifurcation (list of integers).
        """
        return sorted([ei_succ for ei in self.edge_ids_with_paths()
                               for ei_succ in self.edge_successors[ei]
                               if len(self.edge_successors[ei]) > 1])

    @lazy_property
    def edge_predecessors(self):
        """
        Dictionary where edge indices are mapped to a set with the edge predecessors.
        """
        predecessors = {ei: set() for ei in self.edge_ids_with_paths()}
        for pathi, path_data in enumerate(self.rbc_data):
            if len(path_data['time']) >= 2:
                dt = np.median(np.diff(path_data['time']))
                eids = path_data['edgeIndex']
                change_ids = np.where(eids[:-1] != eids[1:])[0]
                for i in change_ids:
                    if path_data['time'][i] + 1.5*dt >= path_data['time'][i+1]:
                        predecessors[eids[i+1]].add(eids[i])
        return predecessors

    @lazy_property
    def edge_successors(self):
        """
        Dictionary where edge indices are mapped to a set with the edge successors.
        """
        successors = {ei: set() for ei in self.edge_ids_with_paths()}
        for pathi, path_data in enumerate(self.rbc_data):
            if len(path_data['time']) >= 2:
                dt = np.median(np.diff(path_data['time']))
                eids = path_data['edgeIndex']
                change_ids = np.where(eids[:-1] != eids[1:])[0]
                for i in change_ids:
                    if path_data['time'][i] + 1.5*dt >= path_data['time'][i+1]:
                        successors[eids[i]].add(eids[i+1])
        return successors

    def edge_ids_with_paths(self, n_path_min=1):
        """
        Return a list of edge indices which have a minimal number of paths crossing them.

        Args:
            n_path_min (int): minimal number of paths required so that an edge be included

        Returns:
            list of int, edge indices
        """
        nPathOnEdges = self.rbc_data.nPathOnEdges
        return sorted([ei for ei in nPathOnEdges if nPathOnEdges[ei] >= n_path_min])

    @lazy_function
    def complete_path_indices_on_edge(self, edge):
        """
        Return the list of complete RBC paths that have a position on the given edge.

        The paths are sorted by the initial time on the given edge.

        The criterion whether a path is complete on an edge is relative to the other paths,
        since it is not possible to rely on edge length (an edge may be partly outside
        the simulation domain). Therefore, an interval is defined based on the averaged
        minimal and maximal s-coordinates of the paths.

        Args:
            edge (int): edge index

        Returns:
            list of int with path indices, sorted by increasing initial time
        """
        path_on_edge_ids = self.rbc_data.edgeIntersectingPathIds(edge)
        min_coords = [min(self.rbc_data.pathCentersOnEdge(pathI, edge))
                      for pathI in path_on_edge_ids]
        max_coords = [max(self.rbc_data.pathCentersOnEdge(pathI, edge))
                      for pathI in path_on_edge_ids]
        smin_bound, smax_bound = self.cropped_s_coord_interval_with_paths(
                                      edge, mean_interval_fraction=0.8)
        filtered_path_ids = [i_path for i_path, smin, smax
                             in zip(path_on_edge_ids, min_coords, max_coords)
                             if smin <= smin_bound and smax >= smax_bound]
        return sorted(filtered_path_ids,
                      key=lambda i: self.rbc_data.pathTimesOnEdge(i, edge)[0])

    def last_complete_path_indices_on_edge(self, edge, n_path):
        """
        Return the list of nPath last complete RBC paths that have a position on the given edge.

        Args:
            edge (int): edge index
            n_path (int): number of indices

        Returns:
            list of int with path indices, sorted by increasing initial time
        """
        path_ids = self.complete_path_indices_on_edge(edge)
        if len(path_ids) < n_path:
            warnings.warn('Only {:d} complete paths on edge {:d}, but {:d} requested'.
                          format(len(path_ids), edge, n_path))
            return path_ids
        return path_ids[-n_path:]

    def complete_path_ids_on_edge_from_time(self, ei, from_time):
        """
        Return a list of RBC path indices that went through the given edge and entered it
        after the given time.

        Args:
            ei (int): edge index
            from_time (float): minimal time for RBCs to enter

        Returns:
            list of int, RBC path indices
        """
        path_ids = self.complete_path_indices_on_edge(ei)
        return [path_i for path_i in path_ids
                if self.rbc_data.pathTimesOnEdge(path_i, ei)[0] >= from_time]

    @lazy_function
    def complete_path_indices(self):
        """
        Return the list of complete RBC paths through the domain

        Returns:
            list of ints
        """
        return [path_i for path_i in range(len(self.rbc_data))
                if all([path_i in self.complete_path_indices_on_edge(ei)
                        for ei in self.rbc_data.edgesOnPath(path_i)])
                   and self.rbc_data.pathEdgeIndices(path_i)[0] in self.inlet_edges()
                   and self.rbc_data.pathEdgeIndices(path_i)[-1] in self.outlet_edges()
                ]

    @lazy_function
    def complete_path_indices_from_time(self, from_time):
        """
        Return the list of complete RBC paths through the domain that entered after the
        given time.

        Args:
            from_time (float): minimal entering time

        Returns:
            list of int, RBC path indices
        """
        path_ids = self.complete_path_indices()
        return [path_i for path_i in path_ids if self.rbc_data.pathTimes(path_i)[0] >= from_time]

    @lazy_function
    def path_indices_on_edge_through_s_coord(self, edge, s):
        """
        Return the list of RBC paths that cross the coordinate s on the given edge.

        The paths are sorted by the initial time on the given edge.

        Args:
            edge (int): edge index
            s (float): s-coordinate

        Returns:
            list of int with path indices, sorted by increasing initial time on the edge
        """
        path_on_edge_ids = self.rbc_data.edgeIntersectingPathIds(edge)
        filtered_path_ids = [i_path for i_path in path_on_edge_ids
                             if len(np.where(self.rbc_data.pathCentersOnEdge(i_path, edge) < s)[0])
                             and len(np.where(self.rbc_data.pathCentersOnEdge(i_path, edge) > s)[0])]
        return sorted(filtered_path_ids,
                      key=lambda i: self.rbc_data.pathTimesOnEdge(i, edge)[0])

    def last_path_indices_on_edge_through_s_coord(self, edge, s, n_path):
        """
        Return the list of the n_path last RBC paths that cross the coordinate s on the given edge.

        Args:
            edge (int): edge index
            s (float): s-coordinate
            n_path (int): number of indices

        Returns:
            list of int with path indices, sorted by increasing initial time
        """
        path_ids = self.path_indices_on_edge_through_s_coord(edge, s)
        if len(path_ids) < n_path:
            warnings.warn('Only {:d} complete paths on edge {:d}, but {:d} requested'.
                          format(len(path_ids), edge, n_path))
            return path_ids
        return path_ids[-n_path:]

    def cropped_s_coord_interval_with_paths(self, ei, mean_interval_fraction=0.8):
        """
        Compute an s-coordinate interval that is covered by all the RBCs that completely
        travel through the given edge in an OpenFOAM simulation.

        Since an edge might be partly outside the OpenFOAM domain, the edge length cannot
        be directly used to compute this interval. Instead, the minimal and maximal
        s-coordinates of all RBCPaths through the given edge are used.

        Args:
            ei (int): edge index
            mean_interval_fraction (float, optional): fraction of the interval formed by the
                averaged minimal and maximal path coordinates

        Returns:
            tuple of floats with lower and upper bound of the interval
        """
        path_on_edge_ids = self.rbc_data.edgeIntersectingPathIds(ei)
        min_coords = [min(self.rbc_data.pathCentersOnEdge(pathI, ei))
                      for pathI in path_on_edge_ids]
        max_coords = [max(self.rbc_data.pathCentersOnEdge(pathI, ei))
                      for pathI in path_on_edge_ids]
        mean_min = float(np.mean(np.array(min_coords)))
        mean_max = float(np.mean(np.array(max_coords)))
        a = 0.5*(1 + mean_interval_fraction)
        b = 0.5*(1 - mean_interval_fraction)
        return a*mean_min + b*mean_max, b*mean_min + a*mean_max

    def s_coord_with_most_paths(self, ei):
        """
        Return the s-coordinate with most RBC paths passing through it on the given edge.

        Args:
            ei (int): edge index

        Returns:
            float, s-coordinate
        """
        path_on_edge_ids = self.rbc_data.edgeIntersectingPathIds(ei)
        min_coords = [min(self.rbc_data.pathCentersOnEdge(pathI, ei))
                      for pathI in path_on_edge_ids]
        max_coords = [max(self.rbc_data.pathCentersOnEdge(pathI, ei))
                      for pathI in path_on_edge_ids]

        # find a common s-coordinate by recursively removing path indices
        def find_common_s_coord(min_values, max_values):
            if max(min_values) <= min(max_values):
                return 0.5*(max(min_values) + min(max_values))
            else:
                ids_to_remove = [np.argmax(min_values), np.argmin(max_values)]
                for i in sorted(ids_to_remove, reverse=True):
                    del min_values[i]
                    del max_values[i]
                return find_common_s_coord(min_values, max_values)
        return find_common_s_coord(min_coords, max_coords)

    def edge_coordinate_to_time(self, s, eI, pathI):
        """
        Return time at which the path with index pathI goes through the given edge coordinate

        Args:
            pathI (int): path index
            eI (int): edge index
            s (float): curvilinear coordinate

        Raises:
            IndexError: if no s-coordinate great than s is found
            ValueError: if no s-coordinate lower than s is found

        Returns:
            float, crossing time of edge coordinate eI, s
        """
        centers = self.rbc_data.pathCentersOnEdge(pathI, eI)
        times = self.rbc_data.pathTimesOnEdge(pathI, eI)
        # if the flow direction is negative, reverse the array order in centers and time
        try:
            if not self.positive_path_flow_on_edge(pathI, eI):
                centers = centers[::-1]
                times = times[::-1]
        except FlowReversalError:
            rbc_name = self.rbc_data.pathIndexToRBCName(pathI)
            print 'Flow reversal for {:s} in edge {:d}'.format(rbc_name, eI)
            print "Using original order of curvilinear coordinates."
            pass
        try:
            i_cross = np.where(centers > s)[0][0]
        except IndexError:
            RBCI = self.rbc_data.pathIndexToRBCName(pathI)
            raise CoordinateNotFoundError("Could not find a coordinate greater than {:g} for {:s} on edge {:d}.".
                                          format(s, RBCI, eI))
        if i_cross == 0:
            RBCI = self.rbc_data.pathIndexToRBCName(pathI)
            raise CoordinateNotFoundError("Could not find a coordinate smaller than {:g} for {:s} on edge {:d}.".
                                          format(s, RBCI, eI))
        i0 = i_cross - 1
        i1 = i_cross
        t0 = times[i0]
        s0 = centers[i0]
        t1 = times[i1]
        s1 = centers[i1]
        return t0 + (t1 - t0)*(s - s0)/(s1 - s0)

    def mean_velocity(self, ei, n_path):
        """
        Return the mean velocity of the last n_path RBC paths on edge ei,
        with positive sign.

        Args:
            ei (int): edge index
            n_path (int): number of paths used for averaging

        Returns:
            float, mean velocity
        """
        s_coord = self.s_coord_with_most_paths(ei)
        path_ids = self.last_path_indices_on_edge_through_s_coord(ei, s_coord, n_path)
        centers = [self.rbc_data.pathCentersOnEdge(path_i, ei) for path_i in path_ids]
        times = [self.rbc_data.pathTimesOnEdge(path_i, ei) for path_i in path_ids]
        rbc_velocities = [(x[-1] - x[0])/(t[-1] - t[0]) for x, t in zip(centers, times)]
        return np.abs(np.mean(rbc_velocities))

    def mean_rbc_flow(self, ei, n_path):
        """
        Return the mean RBC flow of the last n_path RBC paths on edge ei

        Args:
            ei (int): edge index
            n_path (int): number of paths used for averaging

        Returns:
            float, mean velocity
        """
        s_coord = self.s_coord_with_most_paths(ei)
        path_ids = self.last_path_indices_on_edge_through_s_coord(ei, s_coord, n_path)
        mid_times = sorted([self.edge_coordinate_to_time(s_coord, ei, path_i)
                            for path_i in path_ids])
        if len(path_ids) <= 1:
            return 0
        else:
            return (len(path_ids) - 1)/(mid_times[-1] - mid_times[0])


class RBCPathSegmentIndexAdapter(object):
    """
    Class that adapts segment-wise the edge indices of RBC paths.

    Attributes:
        adapter_dict (dict): dictionary with the information on the index adaptation.
            The dict has the following structure:
                eI: [(s_1, segment_id_1), (s_2, segment_id_2), ..., (s_k, segment_id_k)]
            This means that a curvilinear coordinate s with s_i < s < s_{i+1} will be
            assigned to the segment index segment_id_i.
    """

    edge_index_key = 'edgeIndex'
    s_coord_key = 'sCoord'

    def __init__(self, adapter_dict):
        self.adapter_dict = adapter_dict
        self._check_adapter_dict()

    @classmethod
    def from_json(cls, path_to_json):
        fp = open(path_to_json)
        json_data = json.load(fp)
        fp.close()
        adapter_dict = json_data
        adapter_dict = {int(key): value for key, value in adapter_dict.iteritems()}
        return cls(adapter_dict)

    def adapt_edge_indices(self, rbc_data):
        for rbc_fields in rbc_data.rbc_fields.values():
            self.compute_segment_indices_from_graph_coordinate(rbc_fields[self.edge_index_key],
                                                               rbc_fields[self.s_coord_key])
        remove_lazy_properties(rbc_data)

    def segments_on_edge(self, ei):
        """
        Return a list of segment indices that are on the given edge.

        Args:
            ei (int): edge index

        Returns:
            list, segment indices on edge
        """
        if ei in self.adapter_dict:
            return [ei] + [adapter_tuple[1] for adapter_tuple in self.adapter_dict[ei]]

        else:
            return [ei]

    def segment_to_edge_index(self, si):
        """
        Return the edge index where the segment with given index lies.

        Args:
            si (int): segment index

        Returns:
            ei (int), edge index
            index_on_edge,
        """
        return self.segment_to_edge_index_pair(si)[0]

    def segment_to_edge_index_pair(self, si):
        """
        Return a 2-tuple with the edge index and the local segment index on the edge
         where the segment with given index lies.

        Args:
            si (int): segment index

        Returns:
            ei (int), edge index
            index_on_edge (int): local segment index on the edge where it lies
        """
        ei = si
        index_on_edge = 0
        for adapted_ei, adapter_tuples in self.adapter_dict.iteritems():
            for i, (s, segment_i) in enumerate(adapter_tuples):
                if si == segment_i:
                    ei = adapted_ei
                    index_on_edge = i+1
        return ei, index_on_edge

    def compute_segment_indices_from_graph_coordinate(self, ei, scoord):
        """
        Returns the adapted segment index from the edge index and the curvilinear coordinate

        Args:
            ei (np.ndarray): edge indices
            scoord (np.ndarray): curvilinear coordinates

        Returns:
            np.ndarray, array with segment indices
        """
        for adapted_ei, adapter_tuples in self.adapter_dict.iteritems():
            for s, segment_i in adapter_tuples:
                ids = np.where((ei == adapted_ei) & (scoord > s))
                ei[ids] = segment_i
        return ei

    def _check_adapter_dict(self):
        """
        Check the consistency of the adapter dictionary.
        """
        for ei, adapter_tuples in self.adapter_dict.iteritems():
            separators = [s for s, _ in adapter_tuples]
            if not all(separators[i] <= separators[i+1] for i in range(len(separators) - 1)):
                raise ValueError('The limiting s-coordinates are not increasing in edge {:d}'
                                 .format(ei))





