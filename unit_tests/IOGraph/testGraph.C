/*---------------------------------------------------------------------------*\
 
Application
    testGraph

Description
    Test basic features of the Boost Graph Library
    This code is based on the quick tour found in
    http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/quick_tour.html

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace boost;

struct Node
{
    int index;
};

struct CapillarySegment
{
    int index;
    scalar diameter;
    scalar flow;
};

struct CapillaryNetwork
{
    word name;
};

template <class Graph> struct exercise_vertex {
    exercise_vertex(Graph& g_) : g(g_) {}
    //...
    typedef typename graph_traits<Graph>
      ::vertex_descriptor Vertex;
    Graph& g;

    void operator()(const Vertex& v) const
    {
      typedef graph_traits<Graph> GraphTraits;
      typename property_map<Graph, vertex_index_t>::type 
        index = get(vertex_index, g);

      std::cout << "out-edges: ";
      typename GraphTraits::out_edge_iterator out_i, out_end;
      typename GraphTraits::edge_descriptor e;
      for (tie(out_i, out_end) = out_edges(v, g); 
           out_i != out_end; ++out_i) {
        e = *out_i;
        Vertex src = source(e, g), targ = target(e, g);
        std::cout << "(" << index[src] << "," 
                  << index[targ] << "), "
                  << g[e].index  << ", "
                  << g[e].diameter  << ", "
                  << g[e].flow  << "; ";
      }
      std::cout << std::endl;

      // Adjacent vertices without edges
      std::cout << "adjacent vertices: ";
      typename graph_traits<Graph>::adjacency_iterator ai;
      typename graph_traits<Graph>::adjacency_iterator ai_end;
      for (tie(ai, ai_end) = adjacent_vertices(v, g);
           ai != ai_end; ++ai)
        std::cout << index[*ai] <<  " ";
      std::cout << std::endl;

    }
    //...
  };


int main(int argc, char *argv[])
{
    // create a typedef for the Graph type
    typedef adjacency_list
    <
        vecS, vecS, undirectedS, Node, CapillarySegment, CapillaryNetwork
    > Graph;

    // Make convenient labels for the vertices
    enum { A, B, C, D, E, N };
    const int num_vertices = N;
    const char* name = "ABCDE";
    const label vindex_array[] = {234, 12, 50, 30, 59};

    // writing out the edges in the graph
    typedef std::pair<int,int> Edge;
    Edge edge_array[] = 
    { Edge(A,B), Edge(A,D), Edge(C,A), Edge(D,C),
      Edge(C,E), Edge(B,D), Edge(D,E) };
    const int num_edges = sizeof(edge_array)/sizeof(edge_array[0]);
    label eindex_array[]     = { 20, 21, 22, 23, 24, 25, 26};
    scalar diameter_array[] = {6, 5.5, 6, 7.2, 4, 5, 6};
    scalar flow_array[]     = {15.2, 16.3, 12.6, 3.3, 60.3, 1.2, 6.8};

    // declare a graph object
    Graph g(num_vertices);
    g[graph_bundle].name = name;
    
    // add the edges to the graph object
    typedef graph_traits<Graph> ::edge_descriptor EdgeD;
    bool found;
    EdgeD e;
    for (int i = 0; i < num_edges; ++i)
    {
        // also get the edge descriptor    
        tie(e, found) = add_edge(edge_array[i].first, edge_array[i].second, g);
        // add edge properties to the graph
        g[e].index    = eindex_array[i];
        g[e].diameter = diameter_array[i];
        g[e].flow     = flow_array[i];
    }

    // get the property map for vertex indices
    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    IndexMap index = get(vertex_index, g);
    // attach properties to vertices
    typedef graph_traits<Graph>::vertex_iterator vertex_iter;
    std::pair<vertex_iter, vertex_iter> vp;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first)
    {
        g[*vp.first].index = vindex_array[index[*vp.first]];
    }

    // output all vertices
    Info<< "vertices(g) = ";
    for (vp = vertices(g); vp.first != vp.second; ++vp.first)
      std::cout<< index[*vp.first] <<  "(" << g[*vp.first].index
               << ") ";
    std::cout<< std::endl;

    // output all edges
    std::cout << "edges(g) = ";
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
        std::cout << "(" << index[source(*ei, g)] 
                  << "," << index[target(*ei, g)] << ") ";
    std::cout << std::endl;

    // adjacency structure of one vertex
    std::for_each(vertices(g).first, vertices(g).second,
                  exercise_vertex<Graph>(g));


    return 0;
}


