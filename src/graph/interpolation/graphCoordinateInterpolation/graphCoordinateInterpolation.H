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

Class
    Foam::graphCoordinateInterpolation

Description
    Interpolate curvilinear coordinates of a graph given positions at
    different times. Position are given by an edge index and a curvilinear
    coordinate along this edge. The curvilinear coordinate must be between 0
    and the respective edge length.

SourceFiles
    graphCoordinateInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef graphCoordinateInterpolation_H
#define graphCoordinateInterpolation_H

#include "Tuple2.H"
#include "scalarList.H"
#include "labelList.H"

#include "geometricEdgeGraph.H"
#include "graphCoordinateDirected.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;

// Forward declaration of friend functions and operators
class graphCoordinateInterpolation;
Ostream& operator<<(Ostream&, const graphCoordinateInterpolation&);


/*---------------------------------------------------------------------------*\
                         Class graphCoordinateInterpolation Declaration
\*---------------------------------------------------------------------------*/

class graphCoordinateInterpolation
{
    // Private data

        //- Description of data_
        const geometricEdgeGraph& graph_;


    // Private Member Functions
    
        //- Interpolate given ordered segment consisting of two edges.
        //  That is, the edge e1 is (A,B), the edge e2 is (B,C)
        //  The curvilinear coordinates s1 and s2 are given over the whole
        //  segment (A,B,C). 
        //  l1 is the length of e1.
        //  lambda is the interpolation weight.
        graphCoordinateDirected interpolateOrderedSegment
        (
            const label e1,
            const label e2,
            const scalar s1,
            const scalar s2,
            const scalar l1,
            const scalar lambda
        ) const;

        //- Find the weight lambda so that the linearly
        //  interpolated position gcOld + lambda*gcNew is equal
        //  to the given graphCoordinate gc. This method only works when
        //  gcOld and gcNew are on the same edge.
        scalar interpWeightOfCoordinateInSameEdge
        (
            const graphCoordinate& gcTest,
            const graphCoordinate& gcOld,
            const graphCoordinate& gcNew
        ) const;

        scalar interpWeightOfCoordinateUpToNeighbourEdge
        (
            const graphCoordinate& gcTest,
            const graphCoordinate& gcOld,
            const graphCoordinate& gcNew
        ) const;

        //- Disallow default bitwise copy construct
        graphCoordinateInterpolation(const graphCoordinateInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const graphCoordinateInterpolation&);


public:

    // Runtime type information
    ClassName("graphCoordinateInterpolation");

    // Constructors

        //- Construct from components
        graphCoordinateInterpolation(const geometricEdgeGraph&);


    //- Destructor
    ~graphCoordinateInterpolation();


    // Member Functions

        //- Return whether the coordinate gcTest is on the same edge as and
        //  between gc1 and gc2.
        bool isBetweenOnSameEdge
        (
            const graphCoordinate& gcTest,
            const graphCoordinate& gc1,
            const graphCoordinate& gc2
        ) const;

        //- Return whether the coordinate gcTest is between gc1 and gc2, under
        //  the condition that the three points are separated by at most one
        //  vertex.
        bool isBetweenUpToNeighbourEdge
        (
            const graphCoordinate& gcTest,
            const graphCoordinate& gc1,
            const graphCoordinate& gc2
        ) const;

        // Interpolation

            //- Interpolate the graph coordinate at a given time t, as a
            //  function of times, edge indices and curvilinear
            //  coordinates
            graphCoordinateDirected
            interpolate
            (
                const scalarList& times,
                const labelList& edges,
                const scalarList& sCoords,
                const scalar t
            ) const;

            //- Compute the time at which a given graph coordinate is attained,
            //  as a function of times, edge indices and curvilinear
            //  coordinates.
            scalar
            computeGraphCoordinatePassingTime
            (
                const scalarList& times,
                const labelList& edges,
                const scalarList& sCoords,
                const graphCoordinate& gc
            ) const;

    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const graphCoordinateInterpolation&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
