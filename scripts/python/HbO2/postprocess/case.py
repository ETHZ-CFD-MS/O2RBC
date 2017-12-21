"""
Classes for postprocessing of OpenFOAM cases.
"""

from collections import OrderedDict
import copy
import numpy as np
import os
import warnings

from HbO2.postprocess.rbcdata.data import make_rbc_data
from HbO2.postprocess.rbcdata.postprocess import RBCDataPostProcessor, \
                                                 RBCDataGraphPostProcessor, \
                                                 RBCPathSegmentIndexAdapter, \
                                                 FlowReversalError, CoordinateNotFoundError
from HbO2.setup.case import SimulationParametersFactory
from HbO2.setup.utils import isAxisymmetricCase, isGraphCase, isGraphStraightCapillariesCase
from parse.case import time_dirs
from utilities.decorators import lazy_property


class CasePostProcessor(object):
    """
    General post processor for an OpenFOAM case.

    The attributes simParams, rbcData and rbcDataPostProcessor are constructed according to
    the type of OpenFOAM case.

    Results are accessible in string form, based on the class attributes result_output_dict
    which should be redefined by subclasses.

    Attributes:
        case_path (str): path to case to postprocess
        simParams (HbO2Parameters: simulation parameters

    Lazy properties:
        rbc_data (rbc_data): RBC data
        rbcDataPostProcessor (RBCDataPostProcessor)

    Class attributes:
        result_output_dict (OrderedDict): dictionary that defines the postprocessor text output
        output_file_name (str): name of the output file.
    """

    output_file_name = 'results.txt'

    # Dictionary that defines how to extract and display results in the method result_string.
    # The syntax is:
    # dict[headerName] = (methodName, args, formatString), where
    #     headerName (str): name of the header
    #     methodName: (str) method of postprocessor class to call
    #     args (tuple): method arguments
    #     formatString (str): string formatter for output
    result_output_dict = OrderedDict()

    def __init__(self, case_path):
        self.case_path = case_path
        self.simParams = SimulationParametersFactory().make_sim_params(case_path)

    def run(self):
        pass

    def fieldTimeDirs(self):
        """Return the list of time folders where fields have been written."""
        return time_dirs(os.path.join(self.case_path, 'domain'))

    def sValues(self):
        """
        Return the curvilinear coordinates to use on the edge.

        For simulations in graphs, an offset might need to be added to self.xValues.
        """
        if isAxisymmetricCase(self.case_path):
            return self.xValues
        elif isGraphStraightCapillariesCase(self.case_path):
            return self.xValues + self.simParams.sCoordOffset()
        else:
            raise RuntimeError('Case type not supported')

    @lazy_property
    def rbc_data(self):
        return make_rbc_data(self.case_path)

    @lazy_property
    def rbcDataPostProcessor(self):
        if isGraphCase(self.case_path):
            return RBCDataGraphPostProcessor(self.rbc_data)
        elif isAxisymmetricCase(self.case_path):
            return RBCDataPostProcessor(self.rbc_data)

    def result_string(self, result_name, *args, **kwargs):
        """
        Return a string with the a result value.

        The resultName should be a key of the dictionary result_output_dict.

        The method can also be a property of the class. In that case, getattr(self, method)
        does not return a callable object and the value of getattr(self, method) is directly
        used, without using the additional arguments.

        Args:
            result_name (str): name of result to compute
            *args: additional argument to be passed to methods
            **kwargs: dictionary-like arguments to be passed to methods

        Returns:
            string with the result value
        """
        method, firstargs, format_str = self.result_output_dict[result_name]
        if callable(getattr(self, method)):
            result = getattr(self, method)(*(firstargs + args), **kwargs)
        else:
            result = getattr(self, method)
        try:
            return format_str.format(result)
        except ValueError:
            print 'Error with result', result_name, 'while formatting string', \
                  format_str, 'with result', result
            raise

    def output_files_and_result_strings(self, *args, **kwargs):
        """
        Return a dictionary with the output file name as a key and  a list with all strings
        produced by the postprocessor as a value

        Args:
            *args: additional argument to be passed to methods
            **kwargs: dictionary-like arguments to be passed to methods

        Returns:
            dictionary with list of results
        """
        if self.result_output_dict:
            return {self.output_file_name: [self.result_string(result_name, *args, **kwargs)
                                            for result_name in self.result_output_dict.keys()]}
        else:
            return {}

    def output_files_and_result_names(self):
        """
        Return a dictionary with the output file name as a key and  a list with the result names
        as a value

        Returns:
            dictionary with header
        """
        if self.result_output_dict:
            return {self.output_file_name: self.result_output_dict.keys()}
        else:
            return {}


class GraphCasePostProcessor(CasePostProcessor):
    """
    Postprocessor for an axisymmetric simulation.

    Attributes:
        graph_data (GraphCapillaries): graph information
        segment_index_file (str): file name with the segment indices for edges with multiple segments
                                  that cross the domain bounding box
        segment_index_adapter (RBCPathSegmentIndexAdapter): adapter for segment indices
    """

    result_output_dict = OrderedDict()
    segment_index_file = 'segmentIndices.json'

    def __init__(self, case_path):
        super(GraphCasePostProcessor, self).__init__(case_path)
        self.graph_data = self.simParams.graph
        try:
            self.segment_index_adapter = RBCPathSegmentIndexAdapter.from_json \
                (os.path.join(self.case_path, self.segment_index_file))
        except IOError:
            self.segment_index_adapter = RBCPathSegmentIndexAdapter({})
        self.segment_index_adapter.adapt_edge_indices(self.rbc_data)

    def upstream_scoord(self, si):
        return self.coordinate_interval_in_box(si)[0]

    def downstream_scoord(self, si):
        return self.coordinate_interval_in_box(si)[1]

    def scoord_interval_length(self, si):
        ei, local_index = self.segment_index_adapter.segment_to_edge_index_pair(si)
        try:
            bounds = self.coordinate_bounds_in_box(ei)[local_index]
            return abs(bounds[1] - bounds[0])
        except IndexError:
            default_value = 10e-6
            print """Coordinate bounds intersecting bounding box could not be
                     computed for segment {:d}.""".format(si)
            print """Returning the default value {:g}""".format(default_value)
            print "The user should check what is going on."
            return default_value

    def edge_intersects_box(self, si):
        """
        Return whether the edge ei intersects the bounding box.

        Args:
            si (int): segment index

        Returns:
            bool
        """
        return len(self.coordinate_bounds_in_box(si)) > 0

    def coordinate_bounds_in_box(self, si):
        """
        Return the coordinate bounds where edge ei intersects with the domain bounding box.

        An offset can be applied to inflate the bounding box so that an edge segment that
        only slightly leaves the bounding box is not considered to leave the domain.

        Args:
            si (int): segment index

        Returns:
            out (list): list of 2-tuples of floats
        """
        min_point = np.array([self.simParams['xMin'],
                              self.simParams['yMin'],
                              self.simParams['zMin']])
        max_point = np.array([self.simParams['xMax'],
                              self.simParams['yMax'],
                              self.simParams['zMax']])
        ei, local_index = self.segment_index_adapter.segment_to_edge_index_pair(si)
        intervals = self.graph_data.edge_coordinate_intervals_in_box(ei, min_point, max_point,
                                                                     buffer_width=self.graph_data.edge_radius(ei))
        if len(intervals) == 0:
            warnings.warn('No intersection of edge {:d} with the bounding box found'.format(ei))
        if len(intervals) > 1:
            warnings.warn('More than one interval for edge {:d} within the bounding box'.format(ei))
        return intervals

    def coordinate_intervals_in_box(self, si):
        """
        Return a list of ordered coordinate intervals where edge ei intersects with the domain
        bounding box.
        """
        bounds_list = self.coordinate_bounds_in_box(si)
        ei, local_index = self.segment_index_adapter.segment_to_edge_index_pair(si)
        try:
            positive_flow = self.rbcDataPostProcessor.rbc_path_analyzer.positive_flow(ei)
        except FlowReversalError:
            positive_flow = True
        for i, bounds in enumerate(bounds_list):
            if not positive_flow:
                bounds_list[i] = (bounds[1], bounds[0])
        return bounds_list

    def coordinate_interval_in_box(self, si):
        """
        Return the coordinate interval where edge ei intersects with the domain bounding box.
        The coordinate are ordered so that the the first coordinate is upstream and the second
        is downstream.

        Args:
            ei (int): edge index

        Returns:
            2-tuple of floats, ordered
        """
        ei, local_index = self.segment_index_adapter.segment_to_edge_index_pair(si)
        return self.coordinate_intervals_in_box(ei)[local_index]

    def field_on_whole_path(self, field_name, path_i):
        """
        Compute field values on path segments. The values are computed at the locations
        specified by the methods edge_indices_on_path and s_coords_on_path.

        Args:
            field_name (str): name of the field to plot
            path_i (int): index of the RBC path to plot

        Returns:
            list of np.ndarray, arrays with field values on segments
        """
        x_list = self.path_lengths_on_path(path_i)
        s_list = self.s_coords_on_path(path_i)
        eids_list = self.edge_indices_on_path(path_i)
        field_list = []
        for x, s_coords, eids in zip(x_list, s_list, eids_list):
            field = np.zeros(x.shape)
            for ei in set(eids):  # evaluate field in one go on each edge
                ids = np.where(ei == eids)[0]
                field[ids] = self.rbc_data.fieldValueOnEdge(field_name, s_coords[ids], ei, path_i)
            field_list.append(copy.deepcopy(field))
        return field_list

    def path_lengths_on_path(self, path_i):
        """
        Return the list of integrated curvilinear coordinates along all segments of the
        given RBC path

        The first coordinate is set to be zero.
        A list is returned since multiple edge segments through the bounding box can be
        present in one RBC path. The segments are split based on a maximal distance that
        two contiguous points can have on the same edge.

        Args:
            path_i (int): path index

        Returns:
            list of np.ndarray, arrays with integrated curvilinear coordinates
        """
        centers = self.rbc_data.pathCenters(path_i)
        eids = self.rbc_data.pathEdgeIndices(path_i)
        rbc_path_analyzer = self.rbcDataPostProcessor.rbc_path_analyzer
        x_arrays = []
        for index_range in self._path_segment_index_ranges(path_i):
            # interval = self.coordinate_interval_in_box(eids[index_range[0]])
            x = 0.0
            # set the initial coordinate so that the RBC center is at zero when
            # it enters the domain
            # x = centers[index_range[0]] - interval[0]
            x_array = [x]
            for i in index_range[:-1]:
                if eids[i] == eids[i+1]:
                    x += np.abs(centers[i+1] - centers[i])
                else:
                    # add the length contribution before the node
                    if rbc_path_analyzer.positive_path_flow_on_edge(path_i, eids[i]):
                        ei = self.segment_index_adapter.segment_to_edge_index(eids[i])
                        x += self.graph_data['length'][ei] - centers[i]
                    else:
                        x += abs(centers[i])
                    # add the length contribution after the node
                    if rbc_path_analyzer.positive_path_flow_on_edge(path_i, eids[i+1]):
                        x += abs(centers[i+1])
                    else:
                        ei = self.segment_index_adapter.segment_to_edge_index(eids[i+1])
                        x += self.graph_data['length'][ei] - centers[i+1]
                x_array.append(x)
            x_arrays.append(np.array(copy.deepcopy(x_array)))
        return x_arrays

    def times_on_path(self, path_i):
        """
        Return the list of times along all segments of the given RBC path

        Args:
            path_i (int): path index

        Returns:
            list of np.ndarray, arrays with time values
        """
        times = self.rbc_data.pathTimes(path_i)
        return [np.array(times[index_range])
                for index_range in self._path_segment_index_ranges(path_i)]

    def transit_times_on_path(self, path_i):
        transit_times_list = []
        for times, eids in zip(self.times_on_path(path_i), self.edge_indices_on_path(path_i)):
            # interval = self.coordinate_interval_in_box(eids[0])
            # try:
            #     # set the times so that the RBC center enters the domain at time 0
            #     entrance_time = self.rbcDataPostProcessor.rbc_path_analyzer.\
            #         edge_coordinate_to_time(interval[0], eids[0], path_i)
            # except CoordinateNotFoundError:
            #     entrance_time = times[0]
            entrance_time = times[0]
            transit_times_list.append(times - entrance_time)
        return transit_times_list

    def s_coords_on_path(self, path_i):
        """
        Return the list of curvilinear coordinates along all segments of the given RBC path

        Args:
            path_i (int): path index

        Returns:
            list of np.ndarray, arrays with edge curvilinear coordinates
        """
        centers = self.rbc_data.pathCenters(path_i)
        return [np.array(centers[index_range])
                for index_range in self._path_segment_index_ranges(path_i)]

    def edge_indices_on_path(self, path_i):
        """
        Return the list of edge indices along all segments of the given RBC path

        Args:
            path_i (int): path index

        Returns:
            list of np.ndarray, arrays with edge indices
        """
        eids = self.rbc_data.pathEdgeIndices(path_i)
        return [np.array(eids[index_range])
                for index_range in self._path_segment_index_ranges(path_i)]

    def path_segments(self):
        """
        Return a list of 2-tuples (path_i, index_range), where path_i is the path index
        and index_range the range of indices on the current path.

        Returns:
            list of tuples
        """
        segment_list = []
        for path_i in range(len(self.rbc_data)):
            for index_range in self._path_segment_index_ranges(path_i):
                segment_list.append((path_i, index_range))
        return segment_list

    def _path_segment_index_ranges(self, path_i, max_time_diff=1e-2):
        """
        Return a list of ranges that correspond to all path segments on a given path.

        Args:
            path_i (int): path index
            max_time_diff (float): maximal time difference between points on the same segment.

        Returns:
            list of range objects with the index ranges
        """
        times = self.rbc_data.pathTimes(path_i)
        i0 = 0
        ranges = []
        for i in xrange(len(times) - 1):
            if abs(times[i] - times[i+1]) > max_time_diff:
                ranges.append(range(i0, i+1))
                i0 = i+1
        ranges.append(range(i0, len(times)))
        return ranges


class PostProcessorDecorator(CasePostProcessor):
    """
    Decorator for a postprocessor object.
    """

    def __init__(self, decorated):
        if not isinstance(decorated, CasePostProcessor):
            raise TypeError('The constructor argument is not a CasePostProcessor instance.')
        self.decorated = decorated

    def run(self):
        self.decorated.run()

    def output_files_and_result_strings(self, *args, **kwargs):
        inner_dict = self.decorated.output_files_and_result_strings(*args, **kwargs)
        result_dict = inner_dict.copy()
        result_dict.update(super(PostProcessorDecorator, self).output_files_and_result_strings(*args, **kwargs))
        return result_dict

    def output_files_and_result_names(self):
        inner_dict = self.decorated.output_files_and_result_names()
        result_dict = inner_dict.copy()
        result_dict.update(super(PostProcessorDecorator, self).output_files_and_result_names())
        return result_dict

    def __getattr__(self, name):
        return getattr(self.decorated, name)


class GraphPostProcessorDecorator(GraphCasePostProcessor, PostProcessorDecorator):
    """
    Decorator for a graph postprocessor object.
    """

    def __init__(self, decorated):
        if not isinstance(decorated, GraphCasePostProcessor):
            raise TypeError('The constructor argument is not a CasePostProcessor instance.')
        self.decorated = decorated


class PostProcessorWriter(object):
    """
    Writer for a PostProcessorDecorator instance.
    """

    def __init__(self, post_processor):
        self.post_processor = post_processor
        self.separator = '\t'

    def write_results(self, *args, **kwargs):
        files_and_header = self.post_processor.output_files_and_result_names()
        files_and_results = self.post_processor.output_files_and_result_strings(*args, **kwargs)
        file_objects = []
        for file_name in files_and_header.keys():
            file_objects.append(open(file_name, 'w'))

        for fo, header, results in zip(file_objects, files_and_header.values(),
                                       files_and_results.values()):
            fo.write(self.separator.join(header))
            fo.write('\n')
            fo.write(self.separator.join(results))
            fo.write('\n')
            fo.close()



