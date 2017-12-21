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

#include "ellipseAxes.H"

#include "Pair.H"
#include "quaternion.H"
#include "geometricOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ellipseAxes::ellipseAxes()
:
    minorAxis_(vector(VGREAT,VGREAT,VGREAT)),
    majorAxis_(vector(VGREAT,VGREAT,VGREAT)),
    minorRadius_(VGREAT),
    majorRadius_(VGREAT)
{}


Foam::ellipseAxes::ellipseAxes
(
    const vector& minorAxis,
    const vector& majorAxis,
    const scalar minorRadius,
    const scalar majorRadius
)
:
    minorAxis_(minorAxis),
    majorAxis_(majorAxis),
    minorRadius_(minorRadius),
    majorRadius_(majorRadius)
{
    checkAxesAndRadii();
}


Foam::ellipseAxes::ellipseAxes(const ellipseAxes& obj)
:   
    minorAxis_(obj.minorAxis_),
    majorAxis_(obj.majorAxis_),
    minorRadius_(obj.minorRadius_),
    majorRadius_(obj.majorRadius_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ellipseAxes::~ellipseAxes()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
Foam::ellipseAxes::setAxes
(
    const vector& minorAxis, 
    const vector& majorAxis
)
{
    minorAxis_ = minorAxis;
    majorAxis_ = majorAxis;
    checkAxesAndRadii();
}


void
Foam::ellipseAxes::setRadii
(
    const scalar minorRadius, 
    const scalar majorRadius
)
{
    minorRadius_ = minorRadius;
    majorRadius_ = majorRadius;
    checkAxesAndRadii();
}


Foam::ellipseAxes 
Foam::ellipseAxes::interpolate
(
    const scalar lambda, 
    const ellipseAxes& axes2
) const
{
    const ellipseAxes& axes1 = *this;

    // rotation tensor that maps axes1 to axes2
    tensor T = Foam::rotationTensor(Pair<vector>(axes1.minorAxis(), axes1.majorAxis()),
                                    Pair<vector>(axes2.minorAxis(), axes2.majorAxis()));

    if (mag(det(T) + 1) < SMALL)
    {
        WarningIn("ellipseAxes::interpolate(const scalar, const ellipseAxes&)")
            << "The determinant of the rotation tensor is -1, "
            << "something funny may be happening..." << endl;
    }

    // transform the rotation tensor to a quaternion.
    quaternion q(T);

    // use slerp for interpolation
    quaternion interpolatedQ = slerp(quaternion::I, q, lambda);

    // use this quaternion to transform the ellipse axes.
    return ellipseAxes
           (
                interpolatedQ.transform(axes1.minorAxis()),
                interpolatedQ.transform(axes1.majorAxis()),
                (1.0 - lambda)*axes1.minorRadius() + lambda*axes2.minorRadius(),
                (1.0 - lambda)*axes1.majorRadius() + lambda*axes2.majorRadius()
           );
}


void Foam::ellipseAxes::checkAxesAndRadii() const
{
    if (abs(mag(minorAxis_) - 1) > SMALL || abs(mag(majorAxis_) - 1) > SMALL)
    {
        FatalErrorIn("Foam::ellipseAxes::checkAxesAndRadii() const")
            << "The axis vectors are not unit vectors." << nl
            << abort(FatalError);
    }

    if (mag(minorAxis_ & majorAxis_) > SMALL)
    {
        FatalErrorIn("Foam::ellipseAxes::checkAxesAndRadii() const")
            << "The axis vectors are not orthogonal." << nl
            << abort(FatalError);
    }

    if (minorRadius_ > majorRadius_)
    {
        FatalErrorIn("Foam::ellipseAxes::checkAxesAndRadii() const")
            << "The minor radius is larger than the major radius." << nl
            << abort(FatalError);
    }

    if (minorRadius_ <= 0 || majorRadius_ <= 0)
    {
        FatalErrorIn("Foam::ellipseAxes::checkAxesAndRadii() const")
            << "The radii are not positive." << nl
            << abort(FatalError);
    }
}



// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::ellipseAxes::operator=(const ellipseAxes& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::ellipseAxes::operator=(const Foam::ellipseAxes&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    minorAxis_   = rhs.minorAxis();
    majorAxis_   = rhs.majorAxis();
    minorRadius_ = rhs.minorRadius();
    majorRadius_ = rhs.majorRadius();
}

// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream&
Foam::operator<<(Ostream& os, const ellipseAxes& obj)
{
    os << "Minor axis = " << obj.minorAxis() << nl
       << "Major axis = " << obj.majorAxis() << nl
       << "Minor radius = " << obj.minorRadius() << nl
       << "Major radius = " << obj.majorRadius() << endl;

    return os;
}


// ************************************************************************* //
