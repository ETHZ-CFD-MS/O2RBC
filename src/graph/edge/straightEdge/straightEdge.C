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

#include "straightEdge.H"

#include "Pair.H"

#include "polygonalTube.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::straightEdge::straightEdge
(
    const point& start,
    const point& end
)
:
    geometricEdge(start, end)
{}


Foam::straightEdge::straightEdge(const Pair<point>& pp)
:
    geometricEdge(pp)
{}


Foam::straightEdge::straightEdge(Istream& is)
:
    geometricEdge(is)
{}


Foam::straightEdge::straightEdge(const straightEdge& obj)
:
    geometricEdge(obj)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::straightEdge::~straightEdge()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar
Foam::straightEdge::length() const
{
    return mag(endPoint() - startPoint());
}


Foam::point
Foam::straightEdge::pointPosition(const scalar s) const
{
    scalar clampedS = checkAndClampCurvilinearCoord(s);
    scalar lambda = clampedS/length();

    return (1.0 - lambda)*startPoint() + lambda*endPoint();
}

Foam::vector
Foam::straightEdge::tangentVector(const scalar s) const
{
    vector edgeV = endPoint() - startPoint();

    if (edgeV == vector::zero)
    {
        FatalErrorIn
        (
            "Foam::straightEdge::tangentVector(const scalar)"
        ) << "The vertices " << startPoint() << " and " << endPoint()
          << "have the same positions. The edge tangent vector is undefined."
          << nl
          << abort(FatalError);
    }

    return edgeV/mag(edgeV);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::straightEdge::operator=(const straightEdge& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::straightEdge::operator=(const Foam::straightEdge&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    this->geometricEdge::operator=(rhs);


}


Foam::autoPtr<Foam::circularTube>
Foam::straightEdge::makeCircularTube
(
    const scalarList& segmentDiameters,
    const dictionary& tubeOptions
) const
{
    autoPtr<circularTube> pTube
    (
        new polygonalTube
        (
            *this,
            segmentDiameters,
            tubeOptions
        )
    );

    return pTube;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream&
Foam::operator<<(Ostream& os, const straightEdge& obj)
{
    obj.write(os);

    return os;
}


// ************************************************************************* //
