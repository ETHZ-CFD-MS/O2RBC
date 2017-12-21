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
    Foam::straightEdge

Description
    Edge that consists of a straight line.

SourceFiles
    straightEdge.C

\*---------------------------------------------------------------------------*/

#ifndef straightEdge_H
#define straightEdge_H

#include "geometricEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class circularTube;

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class straightEdge;
Istream& operator>>(Istream&, straightEdge&);
Ostream& operator<<(Ostream&, const straightEdge&);


/*---------------------------------------------------------------------------*\
                         Class straightEdge Declaration
\*---------------------------------------------------------------------------*/

class straightEdge
: 
    public geometricEdge
{

public:

    // Constructors

        //- Construct from components
        straightEdge(const point& start, const point& end);

        //- Construct from point pair
        straightEdge(const Pair<point>&);

        //- Construct from Istream
        straightEdge(Istream&);

        //- Construct as copy
        straightEdge(const straightEdge&);


    // Selectors

        //- Select null constructed
        // static autoPtr<straightEdge> New();


    //- Destructor
    virtual ~straightEdge();


    // Member Functions

        // Geometric queries

            virtual scalar length() const;

            //- Return the position of a point on the edge.
            //  The curvilinear coordinate s should be between 0 and the edge
            //  length.
            virtual point pointPosition(const scalar s) const;

            //- Return the unit tangent vector to the edge at a given position.
            //  The curvilinear coordinate s should be between 0 and the edge
            //  length.
            virtual vector tangentVector(const scalar s) const;

            virtual label nSegments() const
            {
                return 2;
            }

        // Factory methods

            // Construct a tube whose type matches the edge type
            virtual autoPtr<circularTube> makeCircularTube
            (
                const scalarList& segmentDiameters,
                const dictionary& tubeOptions
            ) const;

        // Write


    // Member Operators

        void operator=(const straightEdge&);

    // IOstream Operators

        // friend Istream& operator>>(Istream&, straightEdge&);
        friend Ostream& operator<<(Ostream&, const straightEdge&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
