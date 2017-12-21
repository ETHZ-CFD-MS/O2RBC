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
    Foam::geometricEdgeGraph

Description
    Abstract class for graphs in three-dimensional space with geometric
    information about graph edges.

SourceFiles
    geometricEdgeGraph.C

\*---------------------------------------------------------------------------*/

#ifndef geometricEdgeGraph_H
#define geometricEdgeGraph_H

#include "pointGraph.H"

#include "geometricEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class geometricEdgeGraph Declaration
\*---------------------------------------------------------------------------*/

class geometricEdgeGraph
:
    public pointGraph
{
    // Private data
    
        PtrList<geometricEdge> edges_;

    // Private Member Functions
    
        //- Compute the major or "normal" bisection vector
        virtual vector bisectionVector
        (
            const label vertexIndex, 
            const label edgeIndex,
            const bool major
        ) const;

        //- Compute the major or "normal" bisection angle
        virtual scalar bisectionAngle
        (
            const label vertexIndex, 
            const label edgeIndex,
            const bool major
        ) const;

        //- Disallow default bitwise copy construct
        geometricEdgeGraph(const geometricEdgeGraph&);

        //- Disallow default bitwise assignment
        void operator=(const geometricEdgeGraph&);

public:

    // Runtime type information
    ClassName("geometricEdgeGraph");

    // Constructors

        //- Construct from IOobject
        geometricEdgeGraph(const IOobject&);


    //- Destructor
        ~geometricEdgeGraph();


    // Member Functions

        // Access

            //- Constant access to a given edge
            const geometricEdge& edge(const label edgeIndex) const;

        // Geometric queries

            //- Return the unit tangent vector to an edge at a given vertex.
            //  The vector point from the given vertex.
            virtual vector edgeTangentVectorFromVertex
            (
                const label vertexIndex,
                const label edgeIndex
            ) const;

            //- Compute a vertex which is normal to the given edge vector at 
            //  the vertex, and at least to another adjacent edge vector, if any.
            //  If the vertex degree is one, the other vertex of the given edge
            //  is looked up.
            virtual vector normalToVertexPlane
            (
                const label vertexIndex,
                const label edgeIndex
            ) const;

            //- Compute the bisection vector at a given vertex.
            //  When the vertex has degree two, the bisection vector is the vector in
            //  the plane of both edges, such as the angles between itself and the
            //  edge vectors are equal and <= 90 degrees.
            //  When the vertex has degree one, the bisection vector of adjacent
            //  vertices are first looked up. 
            //  When the vertex has degree >= 3, the bisection vector is guaranteed 
            //  to be in the same plane as the edge vector of the given edge.
            vector bisectionVector
            (
                const label vertexIndex, 
                const label edgeIndex
            ) const;

            //- Compute the major bisection vector at a given vertex.
            //  The used vectors are chosen similarly as in bisectionVector(const
            //  label, const label).
            vector majorBisectionVector
            (
                const label vertexIndex, 
                const label edgeIndex
            ) const;

            //- Compute the bisection angle at a given vertex.
            //  The used vectors are chosen similarly as in bisectionVector(const
            //  label, const label).
            scalar bisectionAngle
            (
                const label vertexIndex, 
                const label edgeIndex
            ) const;

            //- Compute the major bisection angle at a given vertex.
            //  The used vectors are chosen similarly as in bisectionVector(const
            //  label, const label).
            scalar majorBisectionAngle
            (
                const label vertexIndex, 
                const label edgeIndex
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
