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

#include "graphCoordinateDirected.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::graphCoordinateDirected::graphCoordinateDirected()
:
    graphCoordinate(),
    goingForward_(true)
{}


Foam::graphCoordinateDirected::graphCoordinateDirected
(
    const label e,
    const scalar sCoord,
    const bool goingForward
)
:
    graphCoordinate(e, sCoord),
    goingForward_(goingForward)
{}


Foam::graphCoordinateDirected::graphCoordinateDirected(const graphCoordinateDirected& gc)
:
    graphCoordinate(gc),
    goingForward_(gc.goingForward_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::graphCoordinateDirected::~graphCoordinateDirected()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::graphCoordinateDirected::operator=(const graphCoordinateDirected& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::graphCoordinateDirected::operator=(const Foam::graphCoordinateDirected&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    (*this).graphCoordinate::operator=(rhs);
    goingForward_ = rhs.goingForward_;
}

// * * * * * * * * * * * * *  IOstream Operators * * * * * * * * * * * * * * //

Istream& Foam::operator>>(Istream& is, graphCoordinateDirected& gc)
{
    is.readBegin("graphCoordinateDirected");
    is >> gc.edgeIndex() >> gc.sCoord() >> gc.goingForward();
    is.readEnd("graphCoordinateDirected");

    is.check("operator>>(Istream&, graphCoordinate&");

    return is;
}

Ostream& Foam::operator<<(Ostream& os, const graphCoordinateDirected& gc)
{
    os  << token::BEGIN_LIST
        << gc.edgeIndex()  << token::SPACE
        << gc.sCoord()     << token::SPACE
        << gc.goingForward()
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
