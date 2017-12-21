import copy
import numpy as np
import os
import warnings

from HbO2.model.kroghSolution import KroghSolution2DCone
from HbO2.postprocess.tissuevolumes import TissueVolumeResults
from HbO2.setup.case import SimulationParametersFactory
from HbO2.setup.graph import directed_integration_ordering
from HbO2.setup.simulationParameters import HbO2ParametersAxisymmetric
from HbO2.setup.utils import GraphInletValue


class HemoglobinGraphIntegrator(object):

    inlet_value_file = 'inletValues.json'

    def __init__(self, case_path, graph, **kwargs):
        """
        Constructor

        Args:
            case_path (str): path to current case
            graph (igraph.Graph): directed graph with flow and geometry information
            **topological_radii (bool): whether to use topological tissue radii (default: True)
            **tissue_volumes_file (str): name of the file with the topological tissue volumes
            **average_hb (bool): whether to average the hemoblogin saturation at converging
                                 bifurcations
        """
        self.case_path = case_path
        self.graph = graph
        self.sim_params = SimulationParametersFactory().make_sim_params(case_path)
        self.inlet_hb = GraphInletValue.from_json(os.path.join(case_path, self.inlet_value_file))
        self.use_topological_radii = kwargs.get('topological_radii', True)
        self.tissue_volumes_results = TissueVolumeResults(case_path)
        self.average_hb = kwargs.get('average_hb', False)
        self.topological_tissue_volumes_dict = \
            {ei: vt for ei, vt in zip(self.tissue_volumes_results.edge_ids(),
                                      self.tissue_volumes_results.topological_tissue_volumes())}
        self.functional_tissue_volumes_dict = \
            {ei: vt for ei, vt in zip(self.tissue_volumes_results.edge_ids(),
                                      self.tissue_volumes_results.functional_tissue_volumes())}
        self.compute()

    def compute(self):
        """
        Compute the hemoglobin saturation at all nodes in the given graph.
        The results are stored in the attribute 'hb' of each vertex and in the attributes
        'upstream_hb' and 'downstream_hb' of each edge.

        This method accomodates distributions of hemoglobin saturation at each node.
        The distribution at a node is described by a list of IntegratedHemoglobinAtVertex objects.
        """
        integration_edge_seq = directed_integration_ordering(self.graph)
        self.graph.es['upstream_hb'] = None
        self.graph.es['downstream_hb'] = None
        self.graph.es['hb_distribution'] = None
        for edge in integration_edge_seq:
            edge['upstream_hb'] = self.upstream_hb(edge)
            self.graph.vs[edge.tuple[0]]['hb'] = edge['upstream_hb']
            if edge['rbc_flow'] == 0.0:
                warnings.warn('The edge {:d} has zero flow'.format(edge['edge_index']))
                edge['rbc_flow'] = 1e-6  # to avoid division by zero
            downstream_hb = []
            for upstream_int_hb in edge['upstream_hb']:
                krogh_sol = KroghSolution2DCone(self._edge_sim_params(edge, upstream_int_hb.hb))
                length = krogh_sol.geometry['domainLength']
                transit_time = length/krogh_sol.vRBC
                hb_v = krogh_sol.saturationAtX(length)
                distal_coord = upstream_int_hb.coordinate + length
                distal_time = upstream_int_hb.transit_time + transit_time
                downstream_hb.append(IntegratedHemoglobinAtVertex(hb_v,
                                                                  weight=upstream_int_hb.weight,
                                                                  coordinate=distal_coord,
                                                                  transit_time=distal_time))
            edge['downstream_hb'] = downstream_hb
            edge['hb_distribution'] = HemoglobinDistributionOnEdge(edge, edge['upstream_hb'],
                                                                   edge['downstream_hb'])
            self.graph.vs[edge.tuple[1]]['hb'] = edge['downstream_hb']
        self._compute_mean_hb()

    def distal_hb(self):
        """
        Return the hemoglobin saturation distribution at each distal node.

        Returns:
            list of list of IntegratedHemoglobinAtVertex objects.
        """
        return self.distal_hb_dict().values()

    def distal_hb_distribution(self, return_edge_indices=False):
        """
        Return distal distribution of hemoglobin saturation value.

        Args:
            return_edge_indices (bool): whether to return the associated edge indices

        Returns: hb, weights, edge_indices (optional)
            hb (np.ndarray): hemoglobin saturation values
            weights (np.ndarray): associated weights (not normalized)
            edge_indices (np.ndarray, optional): associated edge indices
        """
        distal_vs = self.graph.vs(_outdegree=0)
        distal_es = self.graph.es(_target_in=distal_vs.indices)
        hb_values = []
        weights = []
        edge_indices = []
        for e in distal_es:
            hb_distr = e['downstream_hb']
            hb_values.extend([int_hb.hb for int_hb in hb_distr])
            weights.extend([e['rbc_flow']*int_hb.weight for int_hb in hb_distr])
            edge_indices.extend([e['edge_index']]*len(hb_distr))
        if return_edge_indices:
            return np.array(hb_values), np.array(weights), np.array(edge_indices)
        else:
            return np.array(hb_values), np.array(weights)

    def distal_hb_distribution_with_repeat(self):
        """
        Return the distal distribution of hemoglobin saturation values with repeated values
        to account for the flow weighting

        Returns:
            hb (np.ndarray): hemoglobin saturation values, with repetition
        """
        hb_distr, weights = self.distal_hb_distribution(return_edge_indices=False)
        # choose the unit (minimal) weights based as a fraction of the minimal flow
        # weight, but enforce a minimal weight to avoid overflow.
        weight_unit = max(0.1*np.min(weights), 1e-6*np.max(weights))
        hb_repeated = np.zeros((0,))
        for hb, weight in zip(hb_distr, weights):
            n_repeat = np.int(np.round(weight/weight_unit))
            hb_repeated = np.hstack((hb_repeated, [hb]*n_repeat))
        return hb_repeated

    def distal_hb_dict(self):
        """
        Return the hemoglobin saturation value at each distal node.

        Returns:
            dict with edge index as key and a list of IntegratedHemoglobinAtVertex
            objects as value.
        """
        distal_vs = self.graph.vs(_outdegree=0)
        distal_es = self.graph.es(_target_in=distal_vs.indices)
        return {e['edge_index']: self.graph.vs[e.target]['hb'] for e in distal_es}

    def average_distal_hb(self):
        """
        Return the flow-weighted average hemoglobin saturation at the distal nodes

        Returns:
            float, mean distal hemoglobin saturation
        """
        hb, weights = self.distal_hb_distribution()
        return np.average(hb, weights=weights)

    def std_distal_hb(self):
        """
        Return the flow-weighted standard deviation of hemoglobin saturation at the distal nodes

        Returns:
            float
        """
        hb, weights = self.distal_hb_distribution()
        mean = self.average_distal_hb()
        var = np.average((hb - mean)**2, weights=weights)
        return np.sqrt(var)

    def distal_hb_fraction_below(self, threshold):
        hb, weights = self.distal_hb_distribution()
        return sum(weights[hb <= threshold])/sum(weights)

    def distal_rbc_flow(self):
        """
        Return the flow-weighted average hemoglobin saturation at the distal nodes

        Returns:
            float, mean distal hemoglobin saturation
        """
        return np.array(self.distal_rbc_flow_dict().values())

    def distal_rbc_flow_dict(self):
        """
        Return the hemoglobin saturation value at each distal node.

        Returns:
            dict with edge index as key and hemoglobin saturation as value
        """
        distal_vs = self.graph.vs(_outdegree=0)
        distal_es = self.graph.es(_target_in=distal_vs.indices)
        return {e['edge_index']: e['rbc_flow'] for e in distal_es}

    def upstream_hb(self, edge):
        """
        Compute the upstream value of hemoglobin saturation on the given edge

        Args:
            edge (igraph.Edge): edge object

        Returns:
            list of IntegratedHemoglobinAtVertex objects
        """
        v_up = self.graph.vs[edge.tuple[0]]
        if v_up.indegree() == 0:
            hb = self.inlet_hb.inlet_value(edge['edge_index'])
            return [IntegratedHemoglobinAtVertex(hb, weight=1.0, coordinate=0.0,
                                                 transit_time=0.0)]
        else:
            upstream_eids = self.graph.incident(v_up, mode='IN')
            upstream_edges = self.graph.es[upstream_eids]
            if None in upstream_edges['downstream_hb']:
                raise RuntimeError('Upstream value unknown for segment {:d}'.format(edge['edge_index']))
            flow_sum = sum(upstream_edges['rbc_flow'])
            if flow_sum > 0.0:
                if self.average_hb:
                    hb_averaged = sum([e['downstream_hb'][0].hb*e['rbc_flow'] for e in upstream_edges]) \
                                  / sum(upstream_edges['rbc_flow'])
                    coord_averaged = sum([e['downstream_hb'][0].coordinate*e['rbc_flow'] for e in upstream_edges]) \
                                   / sum(upstream_edges['rbc_flow'])
                    time_averaged = sum([e['downstream_hb'][0].transit_time*e['rbc_flow'] for e in upstream_edges]) \
                                  / sum(upstream_edges['rbc_flow'])
                    return [IntegratedHemoglobinAtVertex(hb_averaged, weight=1.0,
                                                         coordinate=coord_averaged,
                                                         transit_time=time_averaged)]
                else:
                    result_list = []
                    for e in upstream_edges:
                        int_hb_list = [IntegratedHemoglobinAtVertex(int_hb.hb,
                                                                    weight=int_hb.weight*e['rbc_flow']/flow_sum,
                                                                    coordinate=int_hb.coordinate,
                                                                    transit_time=int_hb.transit_time)
                                       for int_hb in e['downstream_hb']]
                        result_list.extend(int_hb_list)
                    return result_list
            else:
                result_list = []
                for e in upstream_edges:
                    int_hb_list = [IntegratedHemoglobinAtVertex(0.0,
                                                                weight=int_hb.weight*e['rbc_flow']/flow_sum,
                                                                coordinate=int_hb.coordinate,
                                                                transit_time=int_hb.transit_time)
                                   for int_hb in e['downstream_hb']]
                    result_list.extend(int_hb_list)
                return result_list

    def tissue_radius(self, edge):
        """
        Return the tissue radius for the given edge.

        Args:
            edge (igraph.Edge): edge instance

        Returns:
            float, tissue radius
        """
        if self.use_topological_radii:
            vt = self.topological_tissue_volumes_dict[edge['edge_index']]
            l = edge['length']
            rw = edge['radius_wall']
            return np.sqrt(vt/(l*np.pi) + rw**2) if vt >= 0 else rw
        else:
            i = np.where(self.tissue_volumes_results.edge_ids() == edge['edge_index'])[0][0]
            return self.tissue_volumes_results.functional_tissue_radii_with_negative_values()[i]

    def flow_balance(self):
        """
        Return a dictionary with the flow balance at each connecting node.
        The keys are the vertex indices.
        """
        connecting_vs = self.graph.vs(_degree_gt=1)
        flow_balance_dict = {}
        for v in connecting_vs:
            es_in = self.graph.es[self.graph.incident(v, mode='IN')]
            es_out = self.graph.es[self.graph.incident(v, mode='OUT')]
            print v.index, es_in.indices, es_out.indices
            print es_in['rbc_flow'], es_out['rbc_flow']
            print
            flow_balance_dict[v.index] = sum(es_in['rbc_flow']) - sum(es_out['rbc_flow'])
        return flow_balance_dict

    def _compute_mean_hb(self):
        """
        Compute the averaged hemoglobin saturation at each node and edge from the computed
        distribution.
        """
        def mean_hb(hb_distr):
            return sum([int_hb.hb*int_hb.weight for int_hb in hb_distr]) \
                   /sum([int_hb.weight for int_hb in hb_distr])
        for e in self.graph.es():
            e['upstream_mean_hb'] = mean_hb(e['upstream_hb'])
            e['downstream_mean_hb'] = mean_hb(e['downstream_hb'])
        for v in self.graph.vs():
            v['mean_hb'] = mean_hb(v['hb'])

    def _edge_sim_params(self, edge, hb_up):
        """
        Return a SimulationParameters instance for the given edge

        Args:
            edge (igraph.Edge): edge instance
            hb_up (float): hemoglobin saturation on the upstream side

        Returns:
            HbO2ParametersAxisymmetric instance
        """
        sim_params = copy.deepcopy(self.sim_params)
        sim_params['LDMean'] = edge['ld']
        sim_params['RBCVelocity'] = edge['rbc_velocity']
        sim_params['radiusRBC'] = edge['radius_rbc']
        sim_params['radiusPlasma'] = edge['radius_plasma']
        sim_params['radiusWall'] = edge['radius_wall']
        sim_params['radiusTissueLeft'] = self.tissue_radius(edge)
        sim_params['radiusTissueRight'] = self.tissue_radius(edge)
        sim_params['domainLength'] = edge['length']
        sim_params['HbInlet'] = hb_up
        return HbO2ParametersAxisymmetric(sim_params)


class IntegratedHemoglobinAtVertex(object):

    def __init__(self, hb, weight, coordinate, transit_time):
        self.hb = hb
        self.weight = weight
        self.coordinate = coordinate
        self.transit_time = transit_time


class HemoglobinDistributionOnEdge(object):
    """
    Describes the hemoglobin distribution on a segment, including the weights,
    flow values, coordinate interval and transit times associated to the hemoglobin values.

    Each attribute of the class IntegratedHemoglobinOnSegment is returned in a list
    when accessing the same attribute in this class.
    """

    def __init__(self, edge, upstream_int_hb_list, downstream_int_hb_list):
        self.integrated_hb_list = []
        for up_int_hb, down_int_hb in zip(upstream_int_hb_list, downstream_int_hb_list):
            self.integrated_hb_list.append(
                IntegratedHemoglobinOnSegment(up_int_hb.hb, down_int_hb.hb,
                                              edge['rbc_flow'], up_int_hb.weight,
                                              up_int_hb.coordinate, down_int_hb.coordinate,
                                              up_int_hb.transit_time, down_int_hb.transit_time)
            )

    def __getattr__(self, item):
        try:
            return [getattr(int_hb, item) for int_hb in self.integrated_hb_list]
        except AttributeError:
            raise AttributeError("'{:s}' object has no attribute {:s}".format(type(self), item))

    def __iter__(self):
        return iter(self.integrated_hb_list)

    def __next__(self):
        return next(self.integrated_hb_list)


class IntegratedHemoglobinOnSegment(object):

    def __init__(self, upstream_hb, downstream_hb, total_flow, flow_weight,
                 upstream_path_coord, downstream_path_coord,
                 upstream_transit_time, downstream_transit_time):
        self.upstream_hb = upstream_hb
        self.downstream_hb = downstream_hb
        self.total_flow = total_flow
        self.flow_weight = flow_weight
        self.upstream_path_coord = upstream_path_coord
        self.downstream_path_coord = downstream_path_coord
        self.upstream_transit_time = upstream_transit_time
        self.downstream_transit_time = downstream_transit_time

        self.flow_part = total_flow*flow_weight
