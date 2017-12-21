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
    Foam::graphCoordinateDirected

Description

SourceFiles
    graphCoordinateDirected.C

\*---------------------------------------------------------------------------*/

#ifndef graphCoordinateDirected_H
#define graphCoordinateDirected_H

#include "graphCoordinate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class graphCoordinateDirected;
Istream& operator>>(Istream&, graphCoordinateDirected&);
Ostream& operator<<(Ostream&, const graphCoordinateDirected&);


/*---------------------------------------------------------------------------*\
                         Class graphCoordinateDirected Declaration
\*---------------------------------------------------------------------------*/

class graphCoordinateDirected
:
    public graphCoordinate
{
    // Private data

        //- Whether going forward with respect to curvilinear coordinates
        bool goingForward_;


public:

    // Constructors

        //- Construct null
        graphCoordinateDirected();

        //- Construct from components
        graphCoordinateDirected(const label e, const scalar, const bool forward);

        //- Construct from Istream
        graphCoordinateDirected(Istream&);

        //- Construct as copy
        graphCoordinateDirected(const graphCoordinateDirected&);


    //- Destructor
    ~graphCoordinateDirected();


    // Member Functions

        // Constant access to edge index
        bool goingForward() const
        {
            return goingForward_;
        }

        // Non-constant access to edge index
        bool& goingForward()
        {
            return goingForward_;
        }


    // Member Operators

        void operator=(const graphCoordinateDirected&);


    // IOstream Operators

        friend Istream& operator>>(Istream&, graphCoordinateDirected&);
        friend Ostream& operator<<(Ostream&, const graphCoordinateDirected&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
