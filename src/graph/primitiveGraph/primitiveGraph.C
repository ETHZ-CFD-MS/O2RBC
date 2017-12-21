/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "primitiveGraph.H"

#include <utility>  // for std::pair
#include <iostream> // for std::cout

#include "Time.H"
#include "ListOps.H"


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<typename VProperties, typename EProperties>
typename Foam::primitiveGraph<VProperties, EProperties>::VertexD 
Foam::primitiveGraph<VProperties, EProperties>::findVertex
(
    const label externalVertexIndex
) const
{
    typename boost::unordered_map<label, VertexD>::const_iterator mapEntryP 
                                = vertexIndexToDescriptor_.find(externalVertexIndex);

    if (mapEntryP == vertexIndexToDescriptor_.end())
    {
        FatalErrorIn ("Foam::primitiveGraph::findVertex(const label externalEdgeIndex)")
            << "There is no vertex with index " << externalVertexIndex << "." << nl
            << abort(FatalError);
    }
    return mapEntryP->second;
}


template<typename VProperties, typename EProperties>
typename Foam::primitiveGraph<VProperties, EProperties>::EdgeD
Foam::primitiveGraph<VProperties, EProperties>::findEdge
(
    const label externalEdgeIndex
) const
{
    typename boost::unordered_map<label, EdgeD>::const_iterator mapEntryP 
                                    = edgeIndexToDescriptor_.find(externalEdgeIndex);

    if (mapEntryP == edgeIndexToDescriptor_.end())
    {
        FatalErrorIn ("Foam::primitiveGraph::findEdge(const label externalEdgeIndex)")
            << "There is no edge with index " << externalEdgeIndex << "." << nl
            << abort(FatalError);
    }
    return mapEntryP->second;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<typename VProperties, typename EProperties>
Foam::primitiveGraph<VProperties, EProperties>::primitiveGraph(const IOobject& io)
:
    IOgraph(io),
    vertexIndexToDescriptor_(nVertices()),
    edgeIndexToDescriptor_(nEdges()),
    graph_()
{
    // create vertices and fill map
    forAll(vertexIndices(), i)
    {
        VertexD v = boost::add_vertex(graph_);
        std::pair<label,VertexD> entry(vertexIndices()[i], v);
        vertexIndexToDescriptor_.insert(entry);
    }

    // add the edges to the graph object
    EdgeD e;
    bool found;

    forAll(adjacencyList(), i)
    {
        labelPair listI = adjacencyList()[i];

        // find descriptors of both vertices
        VertexD first  = vertexIndexToDescriptor_[listI.first()];
        VertexD second = vertexIndexToDescriptor_[listI.second()];

        // add edge to graph
        tie(e, found) = boost::add_edge(first, second, graph_);

        // add entry to hash table
        std::pair<label,EdgeD> entry(edgeIndices()[i], e);
        edgeIndexToDescriptor_.insert(entry);
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<typename VProperties, typename EProperties>
Foam::primitiveGraph<VProperties, EProperties>::~primitiveGraph()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<typename VProperties, typename EProperties>
Foam::label 
Foam::primitiveGraph<VProperties, EProperties>::degree
(
    const label externalVertexIndex
) const
{
    VertexD v = findVertex(externalVertexIndex);
    return boost::degree(v, graph_);
}


template<typename VProperties, typename EProperties>
Foam::Pair<typename Foam::primitiveGraph<VProperties, EProperties>::VertexD> 
Foam::primitiveGraph<VProperties, EProperties>::edgeVertices
(
    const label externalEdgeIndex
) const
{
    EdgeD e = findEdge(externalEdgeIndex);
    VertexD src = boost::source(e, graph_);
    VertexD trg = boost::target(e, graph_);

    return Pair<VertexD>(src, trg);
}


template<typename VProperties, typename EProperties>
Foam::List<typename Foam::primitiveGraph<VProperties, EProperties>::VertexD> 
Foam::primitiveGraph<VProperties, EProperties>::adjacentVertices
(
    const label externalVertexIndex
) const
{
    List<VertexD> adjacentList;
    VertexD v = findVertex(externalVertexIndex);
    // add in-vertices
    // inv_adjacency_iterator is not in graph_traits, so just use
    // Graph::inv_adjacency_iterator
    typename Graph::inv_adjacency_iterator iai, iai_end;
    for (tie(iai, iai_end) = inv_adjacent_vertices(v, graph_); iai != iai_end; ++iai)
    {
        adjacentList.append(*iai);
    }
    // add out-vertices
    typename GraphTraits::adjacency_iterator ai, ai_end;
    for (tie(ai, ai_end) = adjacent_vertices(v, graph_); ai != ai_end; ++ai)
    {
        adjacentList.append(*ai);
    }
    return adjacentList;
}

template<typename VProperties, typename EProperties>
Foam::List<typename Foam::primitiveGraph<VProperties, EProperties>::EdgeD> 
Foam::primitiveGraph<VProperties, EProperties>::adjacentEdges
(
    const label externalVertexIndex
) const
{
    List<EdgeD> adjacentList;
    VertexD v = findVertex(externalVertexIndex);
    // add in-edges
    typename Graph::in_edge_iterator in_i, in_end;
    for (tie(in_i, in_end) = in_edges(v, graph_);
         in_i != in_end; ++in_i)
    {
        adjacentList.append(*in_i);
    }
    // add out-edges
    typename GraphTraits::out_edge_iterator out_i, out_end;
    for (tie(out_i, out_end) = out_edges(v, graph_);
         out_i != out_end; ++out_i)
    {
        adjacentList.append(*out_i);
    }
    return adjacentList;
}


template<typename VProperties, typename EProperties>
bool 
Foam::primitiveGraph<VProperties, EProperties>::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<typename VProperties, typename EProperties>
void 
Foam::primitiveGraph<VProperties, EProperties>::operator=(const primitiveGraph& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::primitiveGraph::operator=(const Foam::primitiveGraph&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

template <typename VProperties, typename EProperties>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os, 
    const Foam::primitiveGraph<VProperties, EProperties>& graph
)
{
    graph.write(os);
    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
