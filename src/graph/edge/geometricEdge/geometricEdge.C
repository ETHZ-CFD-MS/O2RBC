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

#include "geometricEdge.H"

#include "Pair.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar
Foam::geometricEdge::checkAndClampCurvilinearCoord
(
    const scalar s
) const
{
    scalar newS = s;
    scalar eLength = length();

    // check that the curvilinear coordinate is positive
    if (s < 0)
    {
        FatalErrorIn("geometricEdge::checkAndClampCurvilinearCoord"
                     "(const scalar)")
            << "The curvilinear coordinate s = " << s << " is negative." << nl
            << abort(FatalError);
    }
    // check that the curvilinear coordinate is less than the edge length.
    // If this happens, its value is clamped.
    if (s > eLength)
    {
        newS = eLength;

        // only issue a warning when the difference is larger than SMALL
        if (s > eLength + SMALL)
        {
            WarningIn("geometricEdge::checkAndClampCurvilinearCoord"
                      "(const scalar)")
                << "The curvilinear coordinate s = " << s
                << " is larger than the edge length (" << eLength << ")." << nl
                << "    Its value was clamped to the edge length." << endl;
        }
    }

    return newS;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geometricEdge::geometricEdge()
:
    startPoint_(point::zero),
    endPoint_(point::zero)
{}


Foam::geometricEdge::geometricEdge
(
    const point& start,
    const point& end
)
:
    startPoint_(start),
    endPoint_(end)
{}

Foam::geometricEdge::geometricEdge(const Pair<point>& pp)
:
    startPoint_(pp.first()),
    endPoint_(pp.second())
{}

Foam::geometricEdge::geometricEdge(Istream& is)
:
    startPoint_(),
    endPoint_()
{
    Pair<point> points(is);
    startPoint_ = points.first();
    endPoint_   = points.second();
}


Foam::geometricEdge::geometricEdge(const geometricEdge& obj)
:
    startPoint_(obj.startPoint_),
    endPoint_(obj.endPoint_)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// Foam::autoPtr<Foam::geometricEdge>
// Foam::geometricEdge::New()
// {
    // return autoPtr<geometricEdge>(new geometricEdge);
// }


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::geometricEdge::~geometricEdge()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::geometricEdge::write(Ostream& os) const
{
    os  << "Geometric edge with start point " << startPoint_
        << " and end point " << endPoint_ << endl;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::geometricEdge::operator=(const geometricEdge& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::geometricEdge::operator=(const Foam::geometricEdge&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    startPoint_ = rhs.startPoint();
    endPoint_   = rhs.endPoint();
}

// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream&
Foam::operator<<(Ostream& os, const geometricEdge& obj)
{
    obj.write(os);
    return os;
}


// ************************************************************************* //
