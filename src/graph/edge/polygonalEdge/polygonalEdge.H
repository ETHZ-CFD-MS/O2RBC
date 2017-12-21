/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::polygonalEdge

Description
    Edge that consist of a polygonal chain.

SourceFiles
    polygonalEdge.C

\*---------------------------------------------------------------------------*/

#ifndef polygonalEdge_H
#define polygonalEdge_H

#include "geometricEdge.H"

#include "List.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class circularTube;

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class polygonalEdge;
Istream& operator>>(Istream&, polygonalEdge&);
Ostream& operator<<(Ostream&, const polygonalEdge&);


/*---------------------------------------------------------------------------*\
                         Class polygonalEdge Declaration
\*---------------------------------------------------------------------------*/

class polygonalEdge
: 
    public geometricEdge
{
    // Private data
    
        // Vertices of the polygonal path (including start and end point)
        List<point> points_;

        // Edge length
        scalar length_;

        // Curvilinear coordinates of the path vertices
        scalarList pathVertexCoords_;

    // Private member functions

        void computeLength();
        void computePathVertexCoords();

public:

    // Constructors

        //- Construct from components
        polygonalEdge(const List<point>&);

        //- Construct from startPoint and endPoint
        polygonalEdge(const point& start, const point& end);

        //- Construct from point pair
        polygonalEdge(const Pair<point>&);

        //- Construct from Istream
        polygonalEdge(Istream&);

        //- Construct as copy
        polygonalEdge(const polygonalEdge&);


    //- Destructor
    virtual ~polygonalEdge();


    // Member Functions

        // Access

            const List<point>& points() const
            {
                return points_;
            }

            label nPathVertices() const
            {
                return points_.size();
            }

            label nSegments() const
            {
                return nPathVertices() - 1;
            }

        // Geometric queries

            virtual scalar length() const
            {
                return length_;
            }

            //- Return the position of a point on the edge.
            virtual point pointPosition(const scalar s) const;

            //- Return the unit tangent vector to the edge at a given position.
            //  The curvilinear coordinate s should be between 0 and the edge
            //  length.
            virtual vector tangentVector(const scalar s) const;

            //- Return the unit tangent vector on a given segment.
            vector tangentVector(const label segmentI) const;

            //- Return the segment end points for a given curvilinear coordinate.
            Pair<point> segmentEndPoints(const scalar s) const;

            //- Return the segment index for a given curvilinear coordinate.
            label segmentIndex(const scalar s) const;

            //- Return the length of the segments of the polygonal path.
            scalarList segmentLengths() const;

            //- Return the curvilinear coordinates of the path vertices.
            const scalarList& pathVertexCoords() const;

            //- Return the unit tangent vector from one of the path vertices.
            vector tangentVectorFromPathVertex
            (
                const label pathVertexI, 
                const bool positiveDirection
            ) const;

        // Factory methods

            // Construct a tube whose type matches the edge type
            virtual autoPtr<circularTube> makeCircularTube
            (
                const scalarList& segmentDiameters,
                const dictionary& tubeOptions
            ) const;

        // Write

            void write(Ostream& os) const;


    // Member Operators

        void operator=(const polygonalEdge&);

    // IOstream Operators

        // friend Istream& operator>>(Istream&, polygonalEdge&);
        friend Ostream& operator<<(Ostream&, const polygonalEdge&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
