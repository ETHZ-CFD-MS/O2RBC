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

#include "geometricEdgeGraph.H"

#include "geometricOps.H"
#include "mathematicalConstants.H"

#include "geometricEdge.H"
#include "polygonalEdge.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(geometricEdgeGraph, 0);
}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

Foam::vector
Foam::geometricEdgeGraph::bisectionVector
(
    const label vertexIndex, 
    const label edgeIndex,
    const bool major
) const
{
    vector v1 = edgeTangentVectorFromVertex(vertexIndex, edgeIndex);
    vector bisect;

    label degree = vertexDegree(vertexIndex);
    if (degree == 1)
    {
        bisect = v1 ^ normalToVertexPlane(vertexIndex, edgeIndex);
    }
    else
    {
        label otherEdgeIndex = otherAdjacentEdgeIndex(vertexIndex, edgeIndex);
        vector v2 = edgeTangentVectorFromVertex(vertexIndex, otherEdgeIndex);
        if (major)
        {
            bisect = Foam::majorBisectionVector(v1, v2);
        }
        else
        {
            bisect = Foam::bisectionVector(v1, v2);
        }
    }

    return bisect/mag(bisect);
}


Foam::scalar
Foam::geometricEdgeGraph::bisectionAngle
(
    const label vertexIndex, 
    const label edgeIndex,
    const bool major
) const
{
    scalar angle = 0;

    label degree = vertexDegree(vertexIndex);

    if (degree == 1)
    {
        angle = 0.5*constant::mathematical::pi;
    }
    else
    {
        vector v1 = edgeTangentVectorFromVertex(vertexIndex, edgeIndex);
        label otherEdgeIndex = otherAdjacentEdgeIndex(vertexIndex, edgeIndex);
        vector v2 = edgeTangentVectorFromVertex(vertexIndex, otherEdgeIndex);

        if (major)
        {
            angle = Foam::majorBisectionAngle(v1, v2);
        }
        else
        {
            angle = Foam::bisectionAngle(v1, v2);
        }

        if (angle < 3.14/4)
        {
            WarningIn("geometricEdgeGraph::bisectionAngle(const label, const label")
                << "The bisection angle at vertex " << vertexIndex
                << " for edge " << edgeIndex << " is small..." << endl;
            angle = 0.25*constant::mathematical::pi;
        }
    }

    return angle;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geometricEdgeGraph::geometricEdgeGraph(const IOobject& io)
:
    pointGraph(io),
    edges_(nEdges())
{
    word edgeType = lookup("edgeType"); 

    if (edgeType == "straight")
    {
        forAll(edges_, eI)
        {
            edges_.set(eI, new polygonalEdge(edgeVertexPositions(eI)));
        }
    }
    else if (edgeType == "polygonal")
    {
        List<List<point> > edgePoints = lookup("edgePoints");
        
        if (edgePoints.size() != nEdges())
        {
            FatalErrorIn("geometricEdgeGraph(const IOobject&)")
                << "The number of elements in edgePoints does not match the "
                << "number of edges." << nl
                << abort(FatalError);
        }

        forAll(edges_, eI)
        {
            edges_.set(eI, new polygonalEdge(edgePoints[eI]));
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::geometricEdgeGraph::~geometricEdgeGraph()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::geometricEdge&
Foam::geometricEdgeGraph::edge(const label edgeIndex) const
{
    return edges_[edgeIndex];
}


Foam::vector
Foam::geometricEdgeGraph::edgeTangentVectorFromVertex
(
    const label vertexIndex,
    const label edgeIndex
) const
{
    checkVertexInEdge(vertexIndex, edgeIndex);

    labelPair vv = edgeVertexIndices(edgeIndex);
    vector result = vector::zero;

    if (vertexIndex == vv.first())
    {
        result = edges_[edgeIndex].tangentVector(0.0);
    }
    else if (vertexIndex == vv.second())
    {
        scalar edgeLength = edges_[edgeIndex].length();
        result = -(edges_[edgeIndex].tangentVector(edgeLength));
    }

    return result;
}


Foam::vector
Foam::geometricEdgeGraph::normalToVertexPlane
(
    const label vertexIndex,
    const label edgeIndex
) const
{
    label degree = vertexDegree(vertexIndex);
    vector normal = vector::zero;

    vector v1 = edgeTangentVectorFromVertex(vertexIndex, edgeIndex);

    if (degree == 1)
    {
        normal = Foam::normalVector(v1);
    }
    else
    {
        label otherEdgeIndex = otherAdjacentEdgeIndex(vertexIndex, edgeIndex);
        vector v2 = edgeTangentVectorFromVertex(vertexIndex, otherEdgeIndex);

        normal = Foam::normalVector(v1, v2);
    }

    return normal/mag(normal);
}


Foam::vector
Foam::geometricEdgeGraph::bisectionVector
(
    const label vertexIndex, 
    const label edgeIndex
) const
{
    return bisectionVector(vertexIndex, edgeIndex, false);
}


Foam::vector
Foam::geometricEdgeGraph::majorBisectionVector
(
    const label vertexIndex, 
    const label edgeIndex
) const
{
    return bisectionVector(vertexIndex, edgeIndex, true);
}


Foam::scalar
Foam::geometricEdgeGraph::bisectionAngle
(
    const label vertexIndex, 
    const label edgeIndex
) const
{
    return bisectionAngle(vertexIndex, edgeIndex, false);
}


Foam::scalar
Foam::geometricEdgeGraph::majorBisectionAngle
(
    const label vertexIndex, 
    const label edgeIndex
) const
{
    return bisectionAngle(vertexIndex, edgeIndex, true);
}


// ************************************************************************* //
