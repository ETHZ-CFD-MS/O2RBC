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

#include "circularTube.H"

#include "Pair.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::circularTube::circularTube
(
    const ellipseAxes& startAxes,
    const ellipseAxes& endAxes
)
:
    startAxes_(startAxes),
    endAxes_(endAxes)
{}


Foam::circularTube::circularTube(const Pair<ellipseAxes>& axesPair)
:
    startAxes_(axesPair.first()),
    endAxes_(axesPair.second())
{}


Foam::circularTube::circularTube(const circularTube& obj)
:
    startAxes_(obj.startAxes_),
    endAxes_(obj.endAxes_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::circularTube::~circularTube()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void 
Foam::circularTube::setStartAxes(const ellipseAxes& startAxes)
{
    startAxes_ = startAxes;
}


void 
Foam::circularTube::setEndAxes(const ellipseAxes& endAxes)
{
    endAxes_ = endAxes;
}


Foam::scalar
Foam::circularTube::meanOuterRadius() const
{
    return 0.5*meanOuterDiameter();
}


Foam::scalar
Foam::circularTube::meanInnerRadius() const
{
    return 0.5*meanInnerDiameter();
}


Foam::scalar
Foam::circularTube::outerRadius(const scalar s) const
{
    return 0.5*outerDiameter(s);
}


Foam::scalar
Foam::circularTube::innerRadius(const scalar s) const
{
    return 0.5*innerDiameter(s);
}


void Foam::circularTube::write(Ostream& os) const
{
    os  << "Circular tube with start axes :" << nl << startAxes_
        << " and end axes " << nl << endAxes_ << endl;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::circularTube::operator=(const circularTube& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::circularTube::operator=(const Foam::circularTube&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    startAxes_ = rhs.startAxes();
    endAxes_   = rhs.endAxes();
}

// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream&
Foam::operator<<(Ostream& os, const circularTube& obj)
{
    obj.write(os);
    return os;
}


// ************************************************************************* //
