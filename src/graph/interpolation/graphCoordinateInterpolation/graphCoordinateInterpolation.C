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

#include "graphCoordinateInterpolation.H"

#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(graphCoordinateInterpolation, 1);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::graphCoordinateDirected
Foam::graphCoordinateInterpolation::interpolateOrderedSegment
(
    const label e1,
    const label e2,
    const scalar s1,
    const scalar s2,
    const scalar l1,
    const scalar lambda
) const
{
    label edgeInterp;
    scalar sInterp;

    // interpolate
    scalar sInterpSegment = (1 - lambda)*s1 + lambda*s2;
    if (sInterpSegment >= l1)
    {
        sInterp = sInterpSegment - l1;
        edgeInterp = e2;
    }
    else
    {
        sInterp = sInterpSegment;
        edgeInterp = e1;
    }
    return graphCoordinateDirected(edgeInterp, sInterp, true);
}


Foam::scalar
Foam::graphCoordinateInterpolation::interpWeightOfCoordinateInSameEdge
(
    const graphCoordinate& gcTest,
    const graphCoordinate& gcOld,
    const graphCoordinate& gcNew
) const
{
    if ((gcOld.sCoord() <= gcTest.sCoord()) && (gcTest.sCoord() <= gcNew.sCoord()))
    {
        return (gcTest.sCoord() - gcOld.sCoord())/(gcNew.sCoord() - gcOld.sCoord());
    }
    else if ((gcNew.sCoord() <= gcTest.sCoord()) && (gcTest.sCoord() <= gcOld.sCoord()))
    {
        return (gcTest.sCoord() - gcNew.sCoord())/(gcOld.sCoord() - gcNew.sCoord());
    }
    else
    {
        FatalErrorIn("graphCoordinateInterpolation::"
                     "interpWeightOfCoordinateInSameEdge("
                     "const graphCoordinate&, const graphCoordinate"
                     "const graphCoordinate&)")
            << "The tested graphCoordinate is not between the provided bounds. "
            << "Please check that the coordinates have been checked using "
            << "isBetweenOnSameEdge(gcTest, gcOld, gcNew) prior to calling "
            << "this method." << nl
            << abort(FatalError);
    }
    return 0.0;
}

Foam::scalar
Foam::graphCoordinateInterpolation::interpWeightOfCoordinateUpToNeighbourEdge
(
    const graphCoordinate& gcTest,
    const graphCoordinate& gcOld,
    const graphCoordinate& gcNew
) const
{
    if (isBetweenOnSameEdge(gcTest, gcOld, gcNew))
    {
        return interpWeightOfCoordinateInSameEdge(gcTest, gcOld, gcNew);
    }
    // else if (gcOld == gcTest)
    // {
        // return 0.0;
    // }
    // else if (gcNew == gcTest)
    // {
        // return 1.0;
    // }
    else
    {
        label eOld  = gcOld.edgeIndex();
        label eNew  = gcNew.edgeIndex();
        label eTest = gcTest.edgeIndex();

        scalar sOld  = gcOld.sCoord();
        scalar sNew  = gcNew.sCoord();
        scalar sTest = gcTest.sCoord();

        // get edge lengths
        scalar lOld = graph_.edge(eOld).length();
        scalar lNew = graph_.edge(eNew).length();

        // find adjacent vertex indices
        labelPair vvOld = graph_.edgeVertexIndices(eOld);
        labelPair vvNew = graph_.edgeVertexIndices(eNew);
        
        // check whether the edges have a common vertex
        if ((vvOld.first()  != vvNew.first()) 
         && (vvOld.first()  != vvNew.second()) 
         && (vvOld.second() != vvNew.first()) 
         && (vvOld.second() != vvNew.second()))
        {
            FatalErrorIn("graphCoordinateInterpolation::interpolate("
                         "const scalarList&, const labelList&, "
                         "const scalarList&, const scalar")
                << "The two edges do not have a vertex in common." << nl
                << abort(FatalError);
        }

        scalar sOldSegment, sNewSegment, sTestSegment;

        // separate cases:
        // let B be the common vertex, A the other vertex of eOld and C the
        // other vertex of eNew.
        // case 1: eOld = (A,B), eNew = (B,C)
        if (vvOld.second() == vvNew.first())
        {
            // find the curvilinear coordinates on the ordered segment (A,B,C)
            sOldSegment = sOld;
            sNewSegment = lOld + sNew;
            if (eTest == eOld)
            {
                sTestSegment = sTest;
            }
            else
            {
                sTestSegment = lOld + sTest;
            }
        }
        // case 2: eOld = (B,A), eNew = (B,C)
        else if (vvOld.first() == vvNew.first())
        {
            // reverse coordinate of first edge
            sOldSegment = lOld - sOld;
            sNewSegment = lOld + sNew;
            if (eTest == eOld)
            {
                sTestSegment = lOld - sTest;
            }
            else
            {
                sTestSegment = lOld + sTest;
            }
        }
        // case 3: eOld = (A,B), eNew = (C,B)
        else if (vvOld.second() == vvNew.second())
        {
            sOldSegment = sOld;
            // reverse coordinate of second edge
            sNewSegment = lOld + lNew - sNew;
            if (eTest == eOld)
            {
                sTestSegment = sTest;
            }
            else
            {
                sTestSegment = lOld + lNew - sTest;
            }
        }
        // case 4: eOld = (B,A), eNew = (C,B)
        else if (vvOld.first() == vvNew.second())
        {
            // reverse coordinate of both edges
            sOldSegment = lOld - sOld;
            sNewSegment = lOld + lNew - sNew;
            if (eTest == eOld)
            {
                sTestSegment = lOld - sTest;
            }
            else
            {
                sTestSegment = lOld + lNew - sTest;
            }
        }
        else
        {
            FatalErrorIn("graphCoordinateInterpolation::"
                         "interpWeightOfCoordinateUpToNeighbourEdge("
                         "const graphCoordinateInterpolation&, "
                         "const graphCoordinateInterpolation&, "
                         "const graphCoordinateInterpolation&)")
                        << "The interpolation weight was not updated, "
                        << "there is certainly a bug in the code." << nl
                        << abort(FatalError);
        }
        return (sTestSegment - sOldSegment)/(sNewSegment - sOldSegment);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::graphCoordinateInterpolation::graphCoordinateInterpolation
(
    const geometricEdgeGraph& graph
)
:
    graph_(graph)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::graphCoordinateInterpolation::~graphCoordinateInterpolation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::graphCoordinateInterpolation::isBetweenOnSameEdge
(
    const graphCoordinate& gcTest,
    const graphCoordinate& gc1,
    const graphCoordinate& gc2
) const
{
    if (gcTest.edgeIndex() == gc1.edgeIndex()
     && gcTest.edgeIndex() == gc2.edgeIndex())
    {
        if (((gc1.sCoord() <= gcTest.sCoord()) && (gcTest.sCoord() <= gc2.sCoord()))
         || ((gc2.sCoord() <= gcTest.sCoord()) && (gcTest.sCoord() <= gc1.sCoord())) )
        {
            return true;
        }
    }

    return false;
}


bool Foam::graphCoordinateInterpolation::isBetweenUpToNeighbourEdge
(
    const graphCoordinate& gcTest,
    const graphCoordinate& gc1,
    const graphCoordinate& gc2
) const
{
    if (isBetweenOnSameEdge(gcTest, gc1, gc2))
    {
        return true;
    }
    else if (gc1 == gcTest || gc2 == gcTest)
    {
        return true;
    }
    else if (gc1.edgeIndex() == gcTest.edgeIndex()
         &&  gcTest.edgeIndex() != gc2.edgeIndex())
    {
        const label commonEdge = gc1.edgeIndex();
        const labelPair vv = graph_.edgeVertexIndices(commonEdge);
        const label vertexOnGcTestSide = 
            (gcTest.sCoord() >= gc1.sCoord() ? vv.second() : vv.first());
        const labelList vEdges = graph_.adjacentEdgeIndices(vertexOnGcTestSide);
        if (findIndex(vEdges, gc2.edgeIndex()) >= 0)
        {
            return true;
        }
    }
    else if (gc2.edgeIndex() == gcTest.edgeIndex()
         &&  gcTest.edgeIndex() != gc1.edgeIndex())
    {
        const label commonEdge = gc2.edgeIndex();
        const labelPair vv = graph_.edgeVertexIndices(commonEdge);
        const label vertexOnGcTestSide = 
            (gcTest.sCoord() >= gc2.sCoord() ? vv.second() : vv.first());
        const labelList vEdges = graph_.adjacentEdgeIndices(vertexOnGcTestSide);
        if (findIndex(vEdges, gc1.edgeIndex()) >= 0)
        {
            return true;
        }
    }
    else if (gc1.edgeIndex() == gc2.edgeIndex()
         &&  gc1.edgeIndex() != gcTest.edgeIndex())
    {
        return false;
    }

    return false;
}


Foam::graphCoordinateDirected
Foam::graphCoordinateInterpolation::interpolate
(
    const scalarList& times,
    const labelList& edges,
    const scalarList& sCoord,
    const scalar t
) const
{
    graphCoordinateDirected gc;

    // find the last index in times which is <= t
    label iOld = -1;
    for(int i=0; i < times.size()-1; ++i)
    {
        if ((times[i] <= t) && (t < times[i+1]))
        {
            iOld = i;
            break;
        }
    }

    if (iOld == -1 || iOld==times.size()-1)
    {
        FatalErrorIn("graphCoordinateInterpolation::interpolate("
                     "const scalarList&, const labelList&, "
                     "const scalarList&, const scalar")
            << "The given time t = " << t << " is outside the range "
            << "of the argument times" << nl
            << abort(FatalError);
    }
    label iNew = iOld+1;

    scalar tOld = times[iOld];
    scalar tNew = times[iNew];

    label eOld = edges[iOld];
    label eNew = edges[iNew];

    scalar sOld = sCoord[iOld];
    scalar sNew = sCoord[iNew];

    // compute the interpolation weight
    scalar lambda = (t - tOld)/(tNew-tOld);

    // if the edges are the same, do a linear interpolation
    if (eOld == eNew)
    {
        scalar sInterp = (1 - lambda)*sOld + lambda*sNew;
        bool forward = (sNew >= sOld);
        gc = graphCoordinateDirected(eOld, sInterp, forward);
    }
    else
    {
        // get edge lengths
        scalar lOld = graph_.edge(eOld).length();
        scalar lNew = graph_.edge(eNew).length();

        // find adjacent vertex indices
        labelPair vvOld = graph_.edgeVertexIndices(eOld);
        labelPair vvNew = graph_.edgeVertexIndices(eNew);
        
        // check whether the edges have a common vertex
        if ((vvOld.first()  != vvNew.first()) 
         && (vvOld.first()  != vvNew.second()) 
         && (vvOld.second() != vvNew.first()) 
         && (vvOld.second() != vvNew.second()))
        {
            FatalErrorIn("graphCoordinateInterpolation::interpolate("
                         "const scalarList&, const labelList&, "
                         "const scalarList&, const scalar")
                << "The two edges do not have a vertex in common." << nl
                << abort(FatalError);
        }

        scalar sOldSegment, sNewSegment;

        // separate cases:
        // let B be the common vertex, A the other vertex of eOld and C the
        // other vertex of eNew.
        // case 1: eOld = (A,B), eNew = (B,C)
        if (vvOld.second() == vvNew.first())
        {
            // find the curvilinear coordinates on the ordered segment (A,B,C)
            sOldSegment = sOld;
            sNewSegment = lOld + sNew;

            // interpolate
            gc = interpolateOrderedSegment(eOld, eNew, sOldSegment,
                                           sNewSegment, lOld, lambda);
        }
        // case 2: eOld = (B,A), eNew = (B,C)
        else if (vvOld.first() == vvNew.first())
        {
            // reverse coordinate of first edge
            sOldSegment = lOld - sOld;
            sNewSegment = lOld + sNew;
            // interpolate
            gc = interpolateOrderedSegment(eOld, eNew, sOldSegment,
                                           sNewSegment, lOld, lambda);
            // if the interpolated point is on the first edge, reverse the
            // coordinate and the direction
            if (gc.edgeIndex() == eOld)
            {
                gc.sCoord() = lOld - gc.sCoord();
                gc.goingForward() = !(gc.goingForward());
            }
        }
        // case 3: eOld = (A,B), eNew = (C,B)
        else if (vvOld.second() == vvNew.second())
        {
            sOldSegment = sOld;
            // reverse coordinate of second edge
            sNewSegment = lOld + lNew - sNew;
            // interpolate
            gc = interpolateOrderedSegment(eOld, eNew, sOldSegment,
                                           sNewSegment, lOld, lambda);
            // if the interpolated point is on the second edge, reverse the
            // coordinate and the direction
            if (gc.edgeIndex() == eNew)
            {
                gc.sCoord() = lNew - gc.sCoord();
                gc.goingForward() = !(gc.goingForward());
            }
        }
        // case 4: eOld = (B,A), eNew = (C,B)
        else if (vvOld.first() == vvNew.second())
        {
            // reverse coordinate of both edges
            sOldSegment = lOld - sOld;
            sNewSegment = lOld + lNew - sNew;
            // interpolate
            gc = interpolateOrderedSegment(eOld, eNew, sOldSegment,
                                           sNewSegment, lOld, lambda);
            // reverse the curvilinear coordinate and the direction
            gc.goingForward() = !(gc.goingForward());
            if (gc.edgeIndex() == eOld)
            {
                gc.sCoord() = lOld - gc.sCoord();
            }
            else
            {
                gc.sCoord() = lNew - gc.sCoord();
            }
        }
        else
        {
            FatalErrorIn("graphCoordinateInterpolation::interpolate("
                         "const scalarList&, const labelList&, "
                         "const scalarList&, const scalar")
                        << "The graph coordinate was not updated, "
                        << "there is certainly a bug in the code." << nl
                        << abort(FatalError);
        }

    }

    return gc;
}

Foam::scalar
Foam::graphCoordinateInterpolation::computeGraphCoordinatePassingTime
(
    const scalarList& times,
    const labelList& edges,
    const scalarList& sCoords,
    const graphCoordinate& gc
) const
{
    scalar passingTime = -VGREAT;
    scalar lambda = 0;
    bool foundCoordinate = false;

    for (int i=0; i < times.size() - 1 && !foundCoordinate; i++)
    {
        graphCoordinate gcStart(edges[i],   sCoords[i]);
        graphCoordinate gcEnd  (edges[i+1], sCoords[i+1]);

        if (isBetweenUpToNeighbourEdge(gc, gcStart, gcEnd))
        {
            lambda = interpWeightOfCoordinateUpToNeighbourEdge(gc, gcStart, gcEnd);
            passingTime = times[i] + lambda*(times[i+1] - times[i]);
            foundCoordinate = true;
        }
    }

    if (!foundCoordinate)
    {
        FatalErrorIn("graphCoordinateInterpolation::computeGraphCoordinatePassingTime("
                     "const scalarList&, const labelList&, "
                     "const scalarList&, const graphCoordinate&")
                    << "The tested graphCoordinate was not found to be between "
                    << "any pair consecutive pair of coordinates." << nl
                    << abort(FatalError);
    }

    return passingTime;
}
    

// ************************************************************************* //
