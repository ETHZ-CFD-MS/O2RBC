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
    Foam::graphCoordinate

Description
    Coordinate for points in a graph. The coordinate consists of the edge index
    and the curvilinear coordinate. The curvilinear coordinate is 0 on the
    first edge vertex and equal to the edge length on the second edge vertex.

SourceFiles
    graphCoordinate.C

\*---------------------------------------------------------------------------*/

#ifndef graphCoordinate_H
#define graphCoordinate_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class graphCoordinate;
Istream& operator>>(Istream&, graphCoordinate&);
Ostream& operator<<(Ostream&, const graphCoordinate&);


/*---------------------------------------------------------------------------*\
                         Class graphCoordinate Declaration
\*---------------------------------------------------------------------------*/

class graphCoordinate
{
    // Private data

        //- Edge index
        label edgeIndex_;

        //- Curvilinear coordinate
        scalar sCoord_;


public:

    // Constructors

        //- Construct null
        graphCoordinate();

        //- Construct from components
        graphCoordinate(const label e, const scalar);

        //- Construct from Istream
        graphCoordinate(Istream&);

        //- Construct as copy
        graphCoordinate(const graphCoordinate&);


    //- Destructor
    ~graphCoordinate();


    // Member Functions

        // Constant access to edge index
        label edgeIndex() const
        {
            return edgeIndex_;
        }

        // Non-constant access to edge index
        label& edgeIndex()
        {
            return edgeIndex_;
        }

        // Constant access to edge index
        scalar sCoord() const
        {
            return sCoord_;
        }

        // Non-constant access to edge index
        scalar& sCoord()
        {
            return sCoord_;
        }


    // Member Operators

        void operator=(const graphCoordinate&);
        bool operator==(const graphCoordinate&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, graphCoordinate&);
        friend Ostream& operator<<(Ostream&, const graphCoordinate&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
