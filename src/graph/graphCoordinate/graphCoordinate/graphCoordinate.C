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

\*---------------------------------------------------------------------------*/

#include "graphCoordinate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::graphCoordinate::graphCoordinate()
:
    edgeIndex_(0),
    sCoord_(0.)
{}


Foam::graphCoordinate::graphCoordinate
(
    const label e,
    const scalar sCoord
)
:
    edgeIndex_(e),
    sCoord_(sCoord)
{}


Foam::graphCoordinate::graphCoordinate(Istream& is)
{
    is >> *this;
}


Foam::graphCoordinate::graphCoordinate(const graphCoordinate& gc)
:
    edgeIndex_(gc.edgeIndex_),
    sCoord_(gc.sCoord_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::graphCoordinate::~graphCoordinate()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

bool Foam::graphCoordinate::operator==(const graphCoordinate& rhs) const
{
    return (edgeIndex() == rhs.edgeIndex() && sCoord() == rhs.sCoord());
}


void Foam::graphCoordinate::operator=(const graphCoordinate& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::graphCoordinate::operator=(const Foam::graphCoordinate&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    edgeIndex_    = rhs.edgeIndex_;
    sCoord_       = rhs.sCoord_;
}

// * * * * * * * * * * * * *  IOstream Operators * * * * * * * * * * * * * * //

Istream& Foam::operator>>(Istream& is, graphCoordinate& gc)
{
    is.readBegin("graphCoordinate");
    is >> gc.edgeIndex_ >> gc.sCoord_;
    is.readEnd("graphCoordinate");

    is.check("operator>>(Istream&, graphCoordinate&");

    return is;
}

Ostream& Foam::operator<<(Ostream& os, const graphCoordinate& gc)
{
    os  << token::BEGIN_LIST
        << gc.edgeIndex_   << token::SPACE
        << gc.sCoord_
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
