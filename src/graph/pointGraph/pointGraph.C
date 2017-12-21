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

#include "pointGraph.H"

#include <utility>  // for std::pair
#include <iostream> // for std::cout

#include "Time.H"
#include "ListOps.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointGraph, 1);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointGraph::pointGraph(const IOobject& io)
:
    primitiveGraph<pointGraphVProp,pointGraphEProp>(io),
    vertexPositions_(lookup("vertexPositions"))
{
    // sanity checks
    if (vertexPositions_.size() != externalVertexIndices().size())
    {
        FatalErrorIn
        (
            "pointGraph::pointGraph(const IOobject& io)"
        )   << "  vertexPositions and vertexIndices do not have the same size" << nl
            << abort(FatalError);
    }

    // add properties to vertices
    forAll(externalVertexIndices(), vI)
    {
        label externalVI = externalVertexIndices()[vI];
        VertexD v = findVertex(externalVI);
        graph_[v].internalIndex = vI;
        graph_[v].externalIndex = externalVI;
        graph_[v].position  = vertexPositions_[vI];
    }

    // add properties to edges
    forAll(externalEdgeIndices(), eI)
    {
        label externalEI = externalEdgeIndices()[eI];
        EdgeD e = findEdge(externalEI);
        graph_[e].internalIndex = eI;
        graph_[e].externalIndex = externalEI;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointGraph::~pointGraph()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::point&
Foam::pointGraph::vertexPosition(const label vertexIndex) const
{
    const label externalVertexIndex = externalVertexIndices()[vertexIndex];
    VertexD v = findVertex(externalVertexIndex);
    return graph_[v].position;
}


Foam::labelPair
Foam::pointGraph::edgeVertexIndices(const label edgeIndex) const
{
    const label externalEdgeIndex = externalEdgeIndices()[edgeIndex];

    Pair<VertexD> vv = edgeVertices(externalEdgeIndex);
    label v1 = graph_[vv.first() ].internalIndex;
    label v2 = graph_[vv.second()].internalIndex;

    return labelPair(v1, v2);
}


Foam::labelList
Foam::pointGraph::adjacentVertexIndices(const label vertexIndex) const
{
    const label externalVertexIndex = externalVertexIndices()[vertexIndex];
    labelList vertexIndices;
    List<primitiveGraph<pointGraphVProp, pointGraphEProp>::VertexD> 
        vertexIndicesD = 
            primitiveGraph<pointGraphVProp, pointGraphEProp>::
                adjacentVertices(externalVertexIndex);

    forAll(vertexIndicesD, i)
    {
        vertexIndices.append(graph_[vertexIndicesD[i]].internalIndex);
    }

    return vertexIndices;
}


Foam::labelList
Foam::pointGraph::adjacentEdgeIndices(const label vertexIndex) const
{
    labelList edgeIndices;

    const label externalVertexIndex = externalVertexIndices()[vertexIndex];
    List<primitiveGraph<pointGraphVProp, pointGraphEProp>::EdgeD> edgeIndicesD
        = primitiveGraph<pointGraphVProp, pointGraphEProp>::
            adjacentEdges(externalVertexIndex);

    forAll(edgeIndicesD, i)
    {
        edgeIndices.append(graph_[edgeIndicesD[i]].internalIndex);
    }

    return edgeIndices;
}


Foam::Pair<Foam::point>
Foam::pointGraph::edgeVertexPositions(const label edgeIndex) const
{
    const label externalEdgeIndex = externalEdgeIndices()[edgeIndex];

    Pair<VertexD> vv = edgeVertices(externalEdgeIndex);
    point v1 = graph_[vv.first() ].position;
    point v2 = graph_[vv.second()].position;

    return Pair<point>(v1, v2);
}


Foam::label
Foam::pointGraph::otherAdjacentEdgeIndex
(
    const label vertexIndex, 
    const label edgeIndex
) const
{
    const labelList& adjacentEdges = adjacentEdgeIndices(vertexIndex);
    if (findIndex(adjacentEdges, edgeIndex) == -1)
    {
        FatalErrorIn
        (
            "Foam::pointGraph::otherAdjacentEdge(const label vertexIndex, const label edgeIndex)"
        )   << "The edge index " << edgeIndex << " is not adjacent to the vertex " 
            << vertexIndex << "." << nl
            << abort(FatalError);
    }

    forAll(adjacentEdges, i)
    {
        label eI = adjacentEdges[i];
        if (eI != edgeIndex)
        {
            return eI;
        }
    }

    return -1;
}

Foam::label
Foam::pointGraph::otherEdgeVertexIndex
(
    const label vertexIndex, 
    const label edgeIndex
) const
{
    labelPair edgeVertices = edgeVertexIndices(edgeIndex);
    if (vertexIndex == edgeVertices.first())
    {
        return edgeVertices.second();
    }
    else if (vertexIndex == edgeVertices.second())
    {
        return edgeVertices.first();
    }
    else
    {
        FatalErrorIn
        (
            "Foam::pointGraph::otherEdgeVertexIndex(const label, const label)"
        )   << "The vertex index " << vertexIndex << " does not belong to edge "
            << edgeIndex << "." << nl
            << abort(FatalError);
    }

    return -1;
}


Foam::label
Foam::pointGraph::vertexDegree(const label vertexIndex) const
{
    label externalVertexIndex = externalVertexIndices()[vertexIndex];
    return primitiveGraph<pointGraphVProp,pointGraphEProp>::
            degree(externalVertexIndex);
}


bool
Foam::pointGraph::isLeafEdge(const label edgeIndex) const
{
    labelPair edgeVertices = edgeVertexIndices(edgeIndex);
    if (vertexDegree(edgeVertices.first() ) == 1 
     || vertexDegree(edgeVertices.second()) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}


Foam::label
Foam::pointGraph::externalToInternalVertexIndex
(
    const label externalVertexIndex
) const
{
    label vertexIndex = findIndex(externalVertexIndices(), externalVertexIndex);

    if (vertexIndex == -1)
    {
        FatalErrorIn
        (
            "Foam::pointGraph::externalToInternalVertexIndex(const label)"
        )   << "The external vertex index " << externalVertexIndex 
            << " does not belong to the graph." << nl
            << abort(FatalError);
    }

    return vertexIndex;
}


Foam::label
Foam::pointGraph::externalToInternalEdgeIndex
(
    const label externalEdgeIndex
) const
{
    label edgeIndex = findIndex(externalEdgeIndices(), externalEdgeIndex);

    if (edgeIndex == -1)
    {
        FatalErrorIn
        (
            "Foam::pointGraph::externalToInternalEdgeIndex(const label)"
        )   << "The external edge index " << externalEdgeIndex 
            << " does not belong to the graph." << nl
            << abort(FatalError);
    }

    return edgeIndex;
}


void
Foam::pointGraph::checkVertexInEdge
(
    const label vertexIndex, 
    const label edgeIndex
) const
{
    labelPair vv = edgeVertexIndices(edgeIndex);
    if (vertexIndex != vv.first() &&
        vertexIndex != vv.second())
    {
        FatalErrorIn("Foam::pointGraph::"
                     "checkVertexInEdge(const label, const label)")
            << "The vertex index " << vertexIndex << " does not belong to "
            << " the edge " << edgeIndex << nl
            << abort(FatalError);
    }
}


bool 
Foam::pointGraph::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const pointGraph& graph)
{
    graph.dictionary::write(os);

    return os;
}


// ************************************************************************* //
