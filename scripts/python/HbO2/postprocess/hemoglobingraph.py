"""
Postprocessing for hemoglobin saturation simulated in a vascular graph using OpenFOAM
"""

from collections import OrderedDict
import copy
import cPickle as pickle
import numpy as np
import os
from scipy.stats.stats import linregress

from HbO2.model.graphintegration import HemoglobinGraphIntegrator
from HbO2.model.tissuevolumes import TissueRadiusFitter
from HbO2.parse.readhemoglobin import TissueVolumesReader
from HbO2.postprocess.case import GraphPostProcessorDecorator
from HbO2.postprocess.rbcdata.postprocess import FlowReversalError
from HbO2.setup.simulationParameters import HbO2ParametersAxisymmetric
from HbO2.setup.utils import GraphInletValue
from utilities.decorators import lazy_function


class HemoglobinPostProcessor(GraphPostProcessorDecorator):
    """
    Postprocessor for hemoglobin saturation in a vascular graph.

    Attributes:
        field_name (str): name of the field to postprocess
        rbc_radius_factor (float): ratio between RBC radius and plasma radius
        wall_radius_factor (float): ratio between endothelium radius and plasma radius
    """

    inlet_value_file = 'inletValues.json'

    def __init__(self, decorated, **kwargs):
        """
        Constructor

        Args:
            decorated (GraphPostProcessor): object to decorate
        """
        super(HemoglobinPostProcessor, self).__init__(decorated)
        self.rbcDataPostProcessor.removePathsWithUndefinedFlowDirection()
        self.field_name = kwargs.get('fieldName', 'Hb_mean')
        self.rbc_radius_factor = kwargs.get('RBCRadiusFactor', 0.8)
        self.rbc_radius_min = kwargs.get('RBCRadiusMin', 1.5e-6)
        self.wall_radius_factor = kwargs.get('wallRadiusFactor', 1.25)
        self.inlet_hb = GraphInletValue.from_json(os.path.join(self.case_path, self.inlet_value_file))


class HemoglobinOnSegmentsPostProcessor(HemoglobinPostProcessor):
    """
    Postprocessor for hemoglobin saturation on capillary segments in a vascular graph.

    Postprocessing is based on edge segments that intersects the bounding box of the
    computational domain, since one edge can enter and leave the domain several times.

    Class attributes:
        output_file_name (str): file name for postprocessing output
        result_output_dict (dict): dictionary with output information
        averaging_modes (list): list with valid averaging modes (strings)

    Attributes:
        averaging_mode (str): averaging mode (must be one of the elements of averaging_modes)
        n_average (int): number of RBC paths used with averaging_mode = 'nRBC'
        averaging_time_start (float): initial time used for averaging with averaging_mode = 'timeStart'
        default_tissue_radius (float): tissue radius used as a first guess for the
                                       computation of functional tissue radii
    """

    output_file_name = 'hemoglobinOnSegmentsResults.txt'
    result_output_dict = OrderedDict()
    result_output_dict['hbMeanUp'] = ('upstream_mean', (), '{:7.5g}')
    result_output_dict['hbMeanDown'] = ('downstream_mean', (), '{:7.5g}')
    result_output_dict['hbMeanDrop'] = ('mean_difference', (), '{:7.5g}')
    result_output_dict['hbMeanSlope'] = ('mean_slope', (), '{:7.5g}')
    result_output_dict['hbStdUp'] = ('upstream_std', (), '{:7.5g}')
    result_output_dict['hbStdDown'] = ('downstream_std', (), '{:7.5g}')
    result_output_dict['hbStdDifference'] = ('std_difference', (), '{:7.5g}')
    result_output_dict['hbStdSlope'] = ('std_slope', (), '{:7.5g}')
    result_output_dict['sInterval'] = ('scoord_interval_length', (), '{:7.5g}')
    result_output_dict['funcTissueRadius'] = ('functional_tissue_radius', (), '{:7.5g}')
    result_output_dict['funcTissueVolume'] = ('functional_tissue_volume', (), '{:7.5g}')
    result_output_dict['funcO2Extraction'] = ('functional_oxygen_extraction_rate', (), '{:7.5g}')
    result_output_dict['topolTissueRadius'] = ('topological_tissue_radius', (), '{:7.5g}')
    result_output_dict['topolTissueVolume'] = ('topological_tissue_volume', (), '{:7.5g}')
    result_output_dict['meanVelocity'] = ('mean_velocity', (), '{:7.5g}')
    result_output_dict['meanRBCFlow'] = ('mean_rbc_flow', (), '{:7.5g}')
    result_output_dict['meanLD'] = ('mean_linear_density', (), '{:7.5g}')
    result_output_dict['meanHT'] = ('mean_hematocrit', (), '{:7.5g}')
    result_output_dict['meanArrivalTT'] = ('mean_arrival_transit_time', (), '{:7.5g}')
    result_output_dict['meanArrivalPath'] = ('mean_arrival_path_length', (), '{:7.5g}')
    result_output_dict['firstTime'] = ('first_used_path_time', (), '{:7.5g}')
    result_output_dict['nPaths'] = ('n_sampled_paths', (), '{:d}')
    result_output_dict['rSquaredLR'] = ('linregress_r_squared', (), '{:7.5g}')
    result_output_dict['slopeLR'] = ('linregress_slope', (), '{:7.5g}')
    result_output_dict['pValueLR'] = ('linregress_p_value', (), '{:7.5g}')

    averaging_modes = ['nRBC', 'timeStart']

    def __init__(self, decorated, **kwargs):
        super(HemoglobinOnSegmentsPostProcessor, self).__init__(decorated, **kwargs)
        self.averaging_mode = kwargs.get('averagingMode', 'nRBC')
        self.n_average = kwargs.get('nAverage', 50)
        self.averaging_time_start = kwargs.get('averagingStart', None)
        self.default_tissue_radius = kwargs.get('defaultTissueRadius', 15e-6)
        if self.averaging_mode not in self.averaging_modes:
            raise ValueError("Invalid averaging mode {:s}".format(self.averaging_mode))
        elif self.averaging_mode == 'timeStart' and self.averaging_time_start == None:
            raise ValueError("Starting time for the averaging unspecified")
        self.tissue_volumes_reader = TissueVolumesReader(self.case_path, 'topologicalTissueVolumes.txt')
        self._check_edge_indices()

    def write_results(self):
        """
        Write the postprocessing results to a file.
        """
        result_header = self.result_output_dict.keys()
        with open(self.output_file_name, 'w') as f:
            f.write('\t'.join(['segment'] + [h for h in result_header]))
            f.write('\n')
            for ei in self.edge_ids():
                print "Postprocessing segment {:d}...".format(ei)
                try:
                    f.write('\t'.join(['{:d}'.format(ei)] +
                                      [self.result_string(resultName, ei)
                                       for resultName in result_header]))
                    f.write('\n')
                except FlowReversalError:
                    print "Skipping segment {:d} due to flow reversal".format(ei)

    @lazy_function
    def edge_ids(self):
        """
        Compute the list of edges used for postprocessing. Only edges that have a
        sufficient number of RBC paths and that intersect the domain bounding box
        are included.

        Returns:
            list with edge indices
        """
        edges_with_paths = self.rbcDataPostProcessor.edgeIdsWithPaths()
        return filter(self.edge_intersects_box, edges_with_paths)

    @lazy_function
    def upstream_mean(self, si):
        if si in self.rbcDataPostProcessor.rbc_path_analyzer.inlet_edges():
            try:
                return self.inlet_hb.inlet_value(si)
            except KeyError:
                pass
        return float(self.rbcDataPostProcessor.fieldAverageOnEdge(
            self.field_name,
            self.upstream_scoord(si), si,
            nAverage=self.n_rbc_average(si)))

    @lazy_function
    def downstream_mean(self, si):
        return float(self.rbcDataPostProcessor.fieldAverageOnEdge(
            self.field_name,
            self.downstream_scoord(si), si,
            nAverage=self.n_rbc_average(si)))

    @lazy_function
    def upstream_std(self, si):
        if si in self.rbcDataPostProcessor.rbc_path_analyzer.inlet_edges():
            try:
                self.inlet_hb.inlet_value(si)
                return 0.0
            except KeyError:
                pass
        return float(self.rbcDataPostProcessor.fieldStdOnEdge(
            self.field_name,
            self.upstream_scoord(si), si,
            nAverage=self.n_rbc_average(si)))

    @lazy_function
    def downstream_std(self, si):
        return float(self.rbcDataPostProcessor.fieldStdOnEdge(
            self.field_name,
            self.downstream_scoord(si), si,
            nAverage=self.n_rbc_average(si)))

    @lazy_function
    def std_difference(self, si):
        """
        Return the drop in standard deviation of hemoglobin saturation
        """
        return self.downstream_std(si) - self.upstream_std(si)

    def hb_difference(self, si):
        """
        Return the hemoglobin drop in segment si for the RBCs used for averaging.
        """
        return self.rbcDataPostProcessor.fieldOnEdge(self.field_name, self.upstream_scoord(si),
                                                     si, self.n_rbc_average(si)) \
             - self.rbcDataPostProcessor.fieldOnEdge(self.field_name, self.downstream_scoord(si),
                                                     si, self.n_rbc_average(si))

    @lazy_function
    def mean_difference(self, si):
        """
        Return the drop in mean hemoglobin saturation.
        """
        return self.downstream_mean(si) - self.upstream_mean(si)

    @lazy_function
    def mean_slope(self, si):
        """
        Return the slope per meter of the mean hemoglobin saturation.
        """
        return (self.downstream_mean(si) - self.upstream_mean(si))/self.scoord_interval_length(si)

    @lazy_function
    def std_slope(self, si):
        """
        Return the slope per meter of the standard deviation of hemoglobin saturation.
        """
        return (self.downstream_std(si) - self.upstream_std(si))/self.scoord_interval_length(si)

    @lazy_function
    def functional_tissue_radius(self, si):
        """
        Compute functional tissue radius based on the mean hemoglobin drop along segment si.

        The functional radius is the radius that a straight Krogh cylinder would have to
        produced the same hemoglobin drop as the one obtained from the simulations.

        If the hemoglobin saturation increases along the edge, the functional tissue radius
        is set to zero.

        Args:
            si (int): segment index

        Returns:
            float
        """
        segment_sim_params = self._segment_sim_params(si)
        fitter = TissueRadiusFitter(segment_sim_params)
        return fitter.fit_mean_tissue_radius(self.upstream_mean(si),
                                             self.downstream_mean(si))

    @lazy_function
    def functional_tissue_volume(self, si):
        """
        Compute functional tissue radius based on the mean hemoglobin drop along segment si.

        The functional radius is the radius that a straight Krogh cylinder would have to
        produced the same hemoglobin drop as the one obtained from the simulations.

        Args:
            si (int): segment index

        Returns:
            float
        """
        segment_sim_params = self._segment_sim_params(si)
        fitter = TissueRadiusFitter(segment_sim_params)
        return fitter.fit_tissue_volume(self.upstream_mean(si),
                                        self.downstream_mean(si))

    @lazy_function
    def functional_oxygen_extraction_rate(self, si):
        segment_sim_params = self._segment_sim_params(si)
        fitter = TissueRadiusFitter(segment_sim_params)
        return fitter.fit_mean_oxygen_extraction_rate(self.upstream_mean(si),
                                                      self.downstream_mean(si))

    def topological_tissue_radius(self, si):
        vol = self.topological_tissue_volume(si)
        length = self.scoord_interval_length(si)
        rw = self.wall_radius_factor*self.graph_data.edge_radius(
                self.segment_index_adapter.segment_to_edge_index(si))
        return np.sqrt(vol/(length*np.pi) + rw**2) if vol >= 0 else 0

    def topological_tissue_volume(self, si):
        i = self.edge_ids().index(si)
        return self.tissue_volumes_reader.topological_tissue_volumes()[i]

    @lazy_function
    def n_sampled_paths(self, si):
        return self.n_rbc_average(si)

    def mean_velocity(self, si):
        return self.rbcDataPostProcessor.rbc_path_analyzer.mean_velocity(si, self.n_rbc_average(si))

    def mean_rbc_flow(self, si):
        return self.rbcDataPostProcessor.rbc_path_analyzer.mean_rbc_flow(si, self.n_rbc_average(si))

    def mean_hematocrit(self, si):
        ei = self.segment_index_adapter.segment_to_edge_index(si)
        vol_rbc = self.simParams['RBCVolume']
        plasma_radius = self.graph_data.edge_radius(ei)
        return vol_rbc*self.mean_rbc_flow(si)/(np.pi*plasma_radius**2*self.mean_velocity(si))

    def mean_linear_density(self, si):
        ei = self.segment_index_adapter.segment_to_edge_index(si)
        vol_rbc = self.simParams['RBCVolume']
        rbc_radius = max(self.rbc_radius_factor*self.graph_data.edge_radius(ei),
                         self.rbc_radius_min)
        l_rbc = vol_rbc/(np.pi*rbc_radius**2)
        return self.mean_rbc_flow(si)*l_rbc/self.mean_velocity(si)

    def mean_arrival_transit_time(self, si):
        path_ids = self.rbcDataPostProcessor.rbc_path_analyzer.\
            last_complete_path_indices_on_edge(si, self.n_rbc_average(si))
        transit_times = np.zeros(len(path_ids),)
        for i, path_i in enumerate(path_ids):
            eids_list = self.edge_indices_on_path(path_i)
            segment_i = [np.any(eids == si) for eids in eids_list].index(True)
            i_on_edge = np.where(eids_list[segment_i] == si)[0][0]
            transit_times[i] = self.transit_times_on_path(path_i)[segment_i][i_on_edge]
        return np.mean(transit_times)

    def mean_arrival_path_length(self, si):
        path_ids = self.rbcDataPostProcessor.rbc_path_analyzer. \
            last_complete_path_indices_on_edge(si, self.n_rbc_average(si))
        arrival_path_lengths = np.zeros(len(path_ids),)
        for i, path_i in enumerate(path_ids):
            eids_list = self.edge_indices_on_path(path_i)
            segment_i = [np.any(eids == si) for eids in eids_list].index(True)
            i_on_edge = np.where(eids_list[segment_i] == si)[0][0]
            arrival_path_lengths[i] = self.path_lengths_on_path(path_i)[segment_i][i_on_edge]
        return np.mean(arrival_path_lengths)

    @lazy_function
    def first_used_path_time(self, si):
        first_path_i = self.rbcDataPostProcessor.rbc_path_analyzer.\
                       last_complete_path_indices_on_edge(si, self.n_rbc_average(si))[0]
        return self.rbcDataPostProcessor.rbc_data.pathTimesOnEdge(first_path_i, si)[0]

    def linregress_hb_drop_with_time_to_previous_rbc(self, si, threshold=np.inf):
        """
        Compute a linear regression of each RBCs hemoglobin saturation drop with the time difference
        to the previous RBC.

        Args:
            si (int): segment index
            threshold (float): maximum value of time difference used for the linear regression

        Returns:
            float tuple, return value of scipy.stats.linregress

        """
        time_difference = self.rbcDataPostProcessor.timeToPreviousRBC(si, self.n_rbc_average(si))
        hb_drop = self.hb_difference(si)[1:]
        filtered_times = time_difference[time_difference < threshold]
        filtered_drops = hb_drop[time_difference < threshold]
        if filtered_times.size:
            return linregress(filtered_times, filtered_drops)
        else:
            return np.nan, np.nan, np.nan, np.nan

    def linregress_slope(self, si):
        return self.linregress_hb_drop_with_time_to_previous_rbc(si)[0]

    def linregress_r_squared(self, si):
        return self.linregress_hb_drop_with_time_to_previous_rbc(si)[2]**2

    def linregress_p_value(self, si):
        return self.linregress_hb_drop_with_time_to_previous_rbc(si)[3]

    def n_rbc_average(self, si):
        """
        Return the number of RBCs used for averaging, depending on the averaging criteria.

        Args:
            si (int): segment index

        Returns:
            int, number of RBCs used for average
        """
        if self.averaging_mode == 'nRBC':
            return self.n_average
        elif self.averaging_mode == 'timeStart':
            return len(self.rbcDataPostProcessor.rbc_path_analyzer.
                       complete_path_ids_on_edge_from_time(si, self.averaging_time_start))

    def _segment_sim_params(self, si):
        """
        Return the simulation parameters with flow and geometric parameters for the given edge.

        Args:
            si (int): segment index

        Returns:
            HbO2ParametersAxisymmetric
        """
        ei = self.segment_index_adapter.segment_to_edge_index(si)
        sim_params = copy.deepcopy(self.simParams)
        sim_params['LDMean'] = self.mean_linear_density(si)
        sim_params['RBCVelocity'] = self.mean_velocity(si)
        sim_params['radiusRBC'] = max(self.rbc_radius_factor*self.graph_data.edge_radius(ei),
                                      self.rbc_radius_min)
        sim_params['radiusPlasma'] = self.graph_data.edge_radius(ei)
        sim_params['radiusWall'] = self.wall_radius_factor*self.graph_data.edge_radius(ei)
        sim_params['radiusTissueLeft'] = self.default_tissue_radius
        sim_params['radiusTissueRight'] = self.default_tissue_radius
        sim_params['domainLength'] = self.scoord_interval_length(si)
        sim_params['HbInlet'] = self.upstream_mean(si)
        return HbO2ParametersAxisymmetric(sim_params)

    def _check_edge_indices(self):
        eids_functional = self.rbcDataPostProcessor.edgeIdsWithPaths()
        eids_topological = self.tissue_volumes_reader.edge_ids()
        if set(eids_functional) < set(eids_topological):
            print 'Edge indices absent from the hemoglobin results: ', \
                  list(set(eids_topological) - set(eids_functional))
        if set(eids_topological) < set(eids_functional):
            print 'Edge indices absent from the topological volumes: ', \
                  list(set(eids_functional) - set(eids_topological))
            raise RuntimeError('Edges from the functional results are absent from the topological volumes')


class HemoglobinOnWholePathsPostProcessor(HemoglobinPostProcessor):
    """
    Postprocessor for hemoglobin saturation along whole RBC paths in a vascular graph.

    TODO: explain how the returned arrays are structured (number of elements)

    Class attributes:
        output_file_name (str): file name for postprocessing output
        result_output_dict (dict): dictionary with output information
        averaging_modes (list): list with valid averaging modes (strings)

    Attributes:
        field_name (str): name of the field to postprocess
    """

    output_file_name = 'hemoglobinOnWholePathsResults.txt'

    result_output_dict = OrderedDict()
    result_output_dict['meanProximalHb'] = ('mean_proximal_hb', (), '{:7.5g}')
    result_output_dict['stdProximalHb'] = ('std_proximal_hb', (), '{:7.5g}')
    result_output_dict['meanDistalHb'] = ('mean_distal_hb', (), '{:7.5g}')
    result_output_dict['stdDistalHb'] = ('std_distal_hb', (), '{:7.5g}')
    result_output_dict['meanHbDrop'] = ('mean_hb_drop', (), '{:7.5g}')
    result_output_dict['stdHbDrop'] = ('std_hb_drop', (), '{:7.5g}')
    result_output_dict['CoVHbDrop'] = ('coeff_variation_hb_drop', (), '{:7.5g}')
    result_output_dict['meanTransitTime'] = ('mean_transit_time', (), '{:7.5g}')
    result_output_dict['stdTransitTime'] = ('std_transit_time', (), '{:7.5g}')
    result_output_dict['CoVTransitTime'] = ('coeff_variation_transit_time', (), '{:7.5g}')
    result_output_dict['nPathSegment'] = ('n_selected_path_segments', (), '{:d}')

    selection_modes = ['nPath', 'firstTime']

    def __init__(self, decorated, **kwargs):
        super(HemoglobinOnWholePathsPostProcessor, self).__init__(decorated, **kwargs)
        self.selection_mode = kwargs.get('selectionMode', 'nPath')
        self.n_path = kwargs.get('nPath', 100)
        self.selection_first_time = kwargs.get('firstTime', None)
        if self.selection_mode not in self.selection_modes:
            raise ValueError("Invalid selection mode {:s}".format(self.selection_mode))

    @lazy_function
    def proximal_hb(self):
        hb_list = []
        for path_i, index_range in self.selected_path_segments():
            hb = self.rbc_data.pathField(path_i, self.field_name)
            hb_list.append(hb[index_range[0]])
        return np.array(hb_list)

    @lazy_function
    def distal_hb(self):
        hb_list = []
        for path_i, index_range in self.selected_path_segments():
            hb = self.rbc_data.pathField(path_i, self.field_name)
            hb_list.append(hb[index_range[-1]])
        return np.array(hb_list)

    def mean_proximal_hb(self):
        return np.mean(self.proximal_hb())

    def std_proximal_hb(self):
        return np.std(self.proximal_hb())

    def mean_distal_hb(self):
        return np.mean(self.distal_hb())

    def std_distal_hb(self):
        return np.std(self.distal_hb())

    def hb_drop(self):
        return self.proximal_hb() - self.distal_hb()

    def mean_hb_drop(self):
        return np.mean(self.hb_drop())

    def std_hb_drop(self):
        return np.std(self.hb_drop())

    def coeff_variation_hb_drop(self):
        return self.std_hb_drop()/self.mean_hb_drop()

    def mean_hb_slope_transit_time(self):
        return self.hb_drop()/self.transit_times()

    def mean_hb_slope_path_length(self):
        return self.hb_drop()/self.transit_path_lengths()

    @lazy_function
    def transit_times(self):
        transit_time_list = []
        for path_i, index_range in self.selected_path_segments():
            times = self.rbc_data.pathTimes(path_i)
            transit_time_list.append(times[index_range[-1]] - times[index_range[0]])
        return np.array(transit_time_list)

    def mean_transit_time(self):
        return np.mean(self.transit_times())

    def std_transit_time(self):
        return np.std(self.transit_times())

    def coeff_variation_transit_time(self):
        return self.std_transit_time()/self.mean_transit_time()

    @lazy_function
    def transit_path_lengths(self):
        transit_path_list = []
        for path_i, index_range in self.selected_path_segments():
            paths = np.hstack(self.path_lengths_on_path(path_i))
            transit_path_list.append(paths[index_range[-1]] - paths[index_range[0]])
        return np.array(transit_path_list)

    def mean_transit_path_length(self):
        return np.mean(self.transit_path_lengths())

    def std_transit_path_length(self):
        return np.std(self.transit_path_lengths())

    def coeff_variation_transit_path_length(self):
        return self.std_transit_path_length()/self.mean_transit_path_length()

    def corr_sim_hb_drop_vs_transit_time(self):
        _, _, rvalue, pvalue, _ = linregress(self.hb_drop(), self.transit_times())
        return rvalue, pvalue

    def corr_sim_hb_drop_vs_transit_path_length(self):
        _, _, rvalue, pvalue, _ = linregress(self.hb_drop(), self.transit_path_lengths())
        return rvalue, pvalue

    def corr_hb_drop_vs_proximal_hb(self):
        _, _, rvalue, pvalue, _ = linregress(self.proximal_hb(), self.hb_drop())
        return rvalue, pvalue

    def corr_mean_hb_slope_transit_time_vs_proximal_hb(self):
        _, _, rvalue, pvalue, _ = linregress(self.proximal_hb(), self.mean_hb_slope_transit_time())
        return rvalue, pvalue

    def corr_mean_hb_slope_path_length_vs_proximal_hb(self):
        _, _, rvalue, pvalue, _ = linregress(self.proximal_hb(), self.mean_hb_slope_path_length())
        return rvalue, pvalue

    def n_selected_path_segments(self):
        return len(self.selected_path_segments())

    def selected_path_segments(self):
        return [path_segment for path_segment in self.path_segments()
                if path_segment[0] in self.selected_path_indices()]

    def selected_path_indices(self):
        if self.selection_mode == 'nPath':
            path_ids = self.rbcDataPostProcessor.rbc_path_analyzer.complete_path_indices()
            if len(path_ids) > self.n_path:
                path_ids = path_ids[:-self.n_path]
            return path_ids
        elif self.selection_mode == 'firstTime':
            return self.rbcDataPostProcessor.rbc_path_analyzer.\
                        complete_path_indices_from_time(self.selection_first_time)


class HemoglobinPostProcessorWithIntegrator(HemoglobinPostProcessor):
    """
    Postprocessor for hemoglobin saturation with an attribute for the integration
    along the network.

    Attributes:
        integrator (HemoglobinGraphIntegrator): integrator for hemoglobin saturation in a graph
    """

    def __init__(self, decorated, **kwargs):
        """
        Constructor
        Args:
            decorated (GraphPostProcessor): object to decorate
            **average_integrated_hb (bool): whether hemoglobin saturation should be averaged
                                            at converging bifurcations
            **integrate_with_topological_radii (bool): whether topological tissue radii should
                                                       be used for integration
            **graph_pickle_name (str): name of the pickle file that contains the graph
        """
        super(HemoglobinPostProcessorWithIntegrator, self).__init__(decorated, **kwargs)
        average_hb = kwargs.get('averageIntegratedHb', False)
        use_topological_radii = kwargs.get('integrateWithTopologicalRadii', True)
        graph = pickle.load(open(kwargs.get('graphPickleName', 'graphForIntegration.pkl'), 'rb'))
        self.integrator = HemoglobinGraphIntegrator(self.case_path, graph,
                                                    average_hb=average_hb,
                                                    topological_radii=use_topological_radii)
