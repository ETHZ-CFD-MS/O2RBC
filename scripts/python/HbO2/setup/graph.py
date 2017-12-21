import igraph
import warnings

from HbO2.postprocess.rbcdata.postprocess import FlowReversalError


def make_graph_from_post_processor(postprocessor):
    """
    Create a directed graph from a postprocessor of an OpenFOAM simulation.

    The edges and vertices hold a dictionary attribute 'index' that corresponds to
    the edges indices on RBC paths. These indices need not be the same as the igraph index.

    Args:
        postprocessor (HemoglobinOnSegmentsPostProcessor): postprocessor object

    Returns:
        igraph.Graph instance
    """
    path_analyzer = postprocessor.rbcDataPostProcessor.rbc_path_analyzer
    graph = igraph.Graph(directed=True)
    for ei in postprocessor.edge_ids():
        try:
            positive_flow = path_analyzer.positive_flow(ei)
        except FlowReversalError:
            warnings.warn("Skipping edge {:d} due to flow reversal".format(ei))
            continue
        # add vertices
        graph.add_vertex()
        graph.add_vertex()
        v0 = graph.vs[graph.vcount()-2]
        v1 = graph.vs[graph.vcount()-1]
        # add oriented edges along flow direction
        oriented_tuple = (v0, v1) if positive_flow else (v1, v0)
        original_ei = postprocessor.segment_index_adapter.segment_to_edge_index(ei)
        graph.add_edge(*oriented_tuple, edge_index=ei, original_edge_index=original_ei)

    # merge the vertices which are connected through RBC paths by a connecting node
    mapping = range(graph.vcount())
    vids_to_delete = set()
    for e in graph.es:
        # build edges that flow in and out of the upstream vertex
        edges_in, edges_out = incident_edges_to_upstream_vertex(e, path_analyzer)
        # update mapping from the list of vertices that build a connecting node
        node_vids = [e_out.tuple[0] for e_out in edges_out] + [e_in.tuple[1] for e_in in edges_in]
        new_node_vi = min(node_vids)  # unambiguous choice of a new vertex id
        for vi in node_vids:
            mapping[vi] = new_node_vi
        vids_to_delete.update(node_vids)
        vids_to_delete.remove(new_node_vi)
    graph.contract_vertices(mapping)
    # if contract_vertices deleted vertices, remove these indices from vids_to_delete
    vids_to_delete.difference_update(range(graph.vcount(), len(mapping)))
    graph.delete_vertices(list(vids_to_delete))

    # add edge attributes
    edge_attribute_names = ['length', 'radius_rbc', 'radius_plasma', 'radius_wall',
                            'ld', 'rbc_velocity', 'rbc_flow']
    for name in edge_attribute_names:
        graph.es[name] = None
    for e in graph.es:
        ei = e['edge_index']
        original_ei = e['original_edge_index']
        e['length'] = postprocessor.scoord_interval_length(ei)
        e['radius_plasma'] = postprocessor.graph_data.edge_radius(original_ei)
        e['radius_rbc'] = max(postprocessor.rbc_radius_factor*e['radius_plasma'],
                              postprocessor.rbc_radius_min)
        e['radius_wall'] = postprocessor.wall_radius_factor*e['radius_plasma']
        e['ld'] = postprocessor.mean_linear_density(ei)
        e['rbc_velocity'] = postprocessor.mean_velocity(ei)
        e['rbc_flow'] = postprocessor.mean_rbc_flow(ei)
    if not graph.is_dag():
        warnings.warn('The produced graph is not a directed acyclic graph', UserWarning)
    return graph


def incident_edges_to_upstream_vertex(e, path_analyzer):
    """
    Compute the edges which are incident to the upstream vertex of the igraph edge
    given as an argument.

    Args:
        e (igraph.Edge): given graph edge
        path_analyzer (RBCPathAnalyzer): RBC path analyzer

    Returns: segments_in, segments_out
        segments_in: set with igraph.Edge instances that go into the upstream vertex
        segments_out: set with igraph.Edge instances that go out of the upstream vertex
    """
    segments_out = {e}
    segments_in = set()
    graph = e.graph
    predecessors = path_analyzer.edge_predecessors[e['edge_index']]
    for ei_pred in predecessors:
        try:
            e_pred = graph.es(edge_index=ei_pred)[0]
            segments_in.add(e_pred)
        except (IndexError, FlowReversalError):
            warnings.warn("Skipping edge predecessor with index {:d}".format(ei_pred))
    for e_in in segments_in:
        for ei_succ in path_analyzer.edge_successors[e_in['edge_index']]:
            try:
                # if multiple segments, select the most upstream one
                e_succ = graph.es(edge_index=ei_succ)[0]
                segments_out.add(e_succ)
            except (IndexError, FlowReversalError):
                warnings.warn("Skipping edge successor with index {:d}".format(ei_succ))
    return segments_in, segments_out


def directed_integration_ordering(graph):
    """
    Compute an edge ordering for integration along the flow direction in a graph.

    The algorithm visits children of a node only if all the incident edges already are
    in the edge sequence. Multi-edges are taken into account: unique elements of node successors
    are visited, since the method successors of igraph returns non-unique elements in case
    of a multi-edge. Additionally, all edges between the visited node and its current successor
    are added to the edge sequence, so every edge in a multi-edge is added.

    Args:
        graph (igraph.Graph): graph object

    Returns:
        list of igraph.Edge objects
    """
    edge_seq = []
    inlet_vs = graph.vs(_indegree=0)
    visited = set()

    def visit_children(node):
        visited.add(node)
        for v in set(node.successors()):
            edge_seq.extend(graph.es(_source=node.index, _target=v.index))
            incident_edges = graph.incident(v.index, mode="IN")
            if set(incident_edges) <= set([e.index for e in edge_seq]):
                visit_children(v)

    for inlet_v in inlet_vs:
        visit_children(inlet_v)

    return edge_seq

def plot_graph(graph):
    layout = graph.layout("kk")
    graph.vs['label'] = graph.vs.indices
    graph.es['label'] = graph.es['edge_index']
    igraph.plot(graph, "graph_custom_indices.png", layout=layout, bbox=(900, 900))
    print 'Wrote graph plot to graph_custom_indices.png.'
    graph.vs['label'] = graph.vs.indices
    graph.es['label'] = graph.es.indices
    igraph.plot(graph, "graph_native_indices.png", layout=layout, bbox=(900, 900))
    print 'Wrote graph plot to graph_native_indices.png.'
