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
    Foam::geometricEdge

Description

SourceFiles
    geometricEdge.C

\*---------------------------------------------------------------------------*/

#ifndef geometricEdge_H
#define geometricEdge_H

#include "point.H"
#include "Pair.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class circularTube;

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class geometricEdge;
Istream& operator>>(Istream&, geometricEdge&);
Ostream& operator<<(Ostream&, const geometricEdge&);


/*---------------------------------------------------------------------------*\
                         Class geometricEdge Declaration
\*---------------------------------------------------------------------------*/

class geometricEdge
{
    // Private data

        //- Start point
        point startPoint_;

        //- End point
        point endPoint_;

protected:

    // Protected member functions

        // Access
        
            // Non-constant access to startPoint.
            point& startPoint()
            {
                return startPoint_;
            }

            // Non-constant access to endPoint.
            point& endPoint()
            {
                return endPoint_;
            }

        scalar checkAndClampCurvilinearCoord(const scalar) const;

public:

    // Constructors

        //- Construct null
        geometricEdge();

        //- Construct from components
        geometricEdge(const point& start, const point& end);

        //- Construct from point pair
        geometricEdge(const Pair<point>&);

        //- Construct from Istream
        geometricEdge(Istream&);

        //- Construct as copy
        geometricEdge(const geometricEdge&);


    // Selectors

        //- Select null constructed
        // static autoPtr<geometricEdge> New();


    //- Destructor
    virtual ~geometricEdge();


    // Member Functions

        // Access
            const point& startPoint() const
            {
                return startPoint_;
            }

            const point& endPoint() const
            {
                return endPoint_;
            }

        // Geometric queries

            virtual scalar length() const = 0;

            //- Return the position of a point on the edge.
            //  The curvilinear coordinate s should be between 0 and the edge
            //  length.
            virtual point pointPosition(const scalar s) const = 0;

            //- Return the unit tangent vector to the edge at a given position.
            //  The curvilinear coordinate s should be between 0 and the edge
            //  length.
            virtual vector tangentVector(const scalar s) const = 0;

            //- Return the number of segments
            virtual label nSegments() const = 0;

        // Factory methods

            // Construct a tube whose type matches the edge type
            virtual autoPtr<circularTube> makeCircularTube
            (
                const scalarList& segmentDiameters,
                const dictionary& tubeOptions
            ) const = 0;

        // Write

            virtual void write(Ostream&) const;


    // Member Operators

        void operator=(const geometricEdge&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        // friend Istream& operator>>(Istream&, geometricEdge&);
        friend Ostream& operator<<(Ostream&, const geometricEdge&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
