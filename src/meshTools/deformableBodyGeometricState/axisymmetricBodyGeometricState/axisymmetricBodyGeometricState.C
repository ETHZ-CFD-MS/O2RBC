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

#include "axisymmetricBodyGeometricState.H"
#include "addToRunTimeSelectionTable.H"

#include "dictionaryEntry.H"
#include "pointField.H"
#include "transform.H"
#include "transformField.H"
#include "volFields.H"
#include "quaternion.H"
#include "septernion.H"

#include "geometricOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(axisymmetricBodyGeometricState, 0);

    addToRunTimeSelectionTable
    (
        deformableBodyGeometricState,
        axisymmetricBodyGeometricState,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        deformableBodyGeometricState,
        axisymmetricBodyGeometricState,
        mesh
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::axisymmetricBodyGeometricState::axisymmetricBodyGeometricState
(
    const point& center,
    const vector& axis,
    const scalar diameter,
    const scalar length
)
:
    deformableBodyGeometricState(center),
    axis_(axis),
    diameter_(diameter),
    length_(length)
{}


Foam::axisymmetricBodyGeometricState::axisymmetricBodyGeometricState
(
    const dictionary& dict
)
:
    deformableBodyGeometricState(dict),
    axis_(dict.lookup("axis")),
    diameter_(readScalar(dict.lookup("diameter"))),
    length_(readScalar(dict.lookup("length")))
{}


Foam::axisymmetricBodyGeometricState::axisymmetricBodyGeometricState
(
    const dictionary& dict,
    const polyMesh& mesh
)
:
    deformableBodyGeometricState(dict),
    axis_(dict.lookup("axis")),
    diameter_(mesh.bounds().span().y()),
    length_(mesh.bounds().span().x())
{
    if (dict.found("diameter") || dict.found("length"))
    {
        WarningIn
        (
            "axisymmetricBodyGeometricState::"
            "axisymmetricBodyGeometricState(const dictionary&, "
            "const polyMesh& mesh"
        ) << "The dictionary entries for the diameter and the length "
          << "were ignored and obtained from the mesh "
          << mesh.name() << " instead." << endl;
    }
}


Foam::axisymmetricBodyGeometricState::axisymmetricBodyGeometricState
(
    Istream& is
)
:
    deformableBodyGeometricState(is),
    axis_(is),
    diameter_(readScalar(is)),
    length_(readScalar(is))
{}


Foam::axisymmetricBodyGeometricState::axisymmetricBodyGeometricState
(
    const axisymmetricBodyGeometricState& obj
)
:
    deformableBodyGeometricState(obj),
    axis_(obj.axis_),
    diameter_(obj.diameter_),
    length_(obj.length_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::boundBox
Foam::axisymmetricBodyGeometricState::bounds() const
{
    // compute the angles between the body axis and the xy, xz and yz-planes
    scalar lengthXY = mag(vector(axis_.x(), axis_.y(),       0.0));
    scalar lengthXZ = mag(vector(axis_.x(),       0.0, axis_.z()));
    scalar lengthYZ = mag(vector(      0.0, axis_.y(), axis_.z()));
    scalar alphaXY  = atan2(axis_.z(), lengthXY);
    scalar alphaXZ  = atan2(axis_.y(), lengthXZ);
    scalar alphaYZ  = atan2(axis_.x(), lengthYZ);

    // compute the positions of the axis end points
    point vertexA = center() + 0.5*length_*axis_/mag(axis_);
    point vertexB = center() - 0.5*length_*axis_/mag(axis_);

    scalar right  = std::max(vertexA.x() + 0.5*diameter_*cos(alphaYZ),
                             vertexB.x() + 0.5*diameter_*cos(alphaYZ));
    scalar left   = std::min(vertexA.x() - 0.5*diameter_*cos(alphaYZ),
                             vertexB.x() - 0.5*diameter_*cos(alphaYZ));
    scalar top    = std::max(vertexA.y() + 0.5*diameter_*cos(alphaXZ),
                             vertexB.y() + 0.5*diameter_*cos(alphaXZ));
    scalar bottom = std::min(vertexA.y() - 0.5*diameter_*cos(alphaXZ),
                             vertexB.y() - 0.5*diameter_*cos(alphaXZ));
    scalar front  = std::max(vertexA.z() + 0.5*diameter_*cos(alphaXY),
                             vertexB.z() + 0.5*diameter_*cos(alphaXY));
    scalar back   = std::min(vertexA.z() - 0.5*diameter_*cos(alphaXY),
                             vertexB.z() - 0.5*diameter_*cos(alphaXY));

    return boundBox(point(left, bottom, back), point(right, top, front));
}


Foam::pointField 
Foam::axisymmetricBodyGeometricState::transformPoints
(
    const pointField& undisplacedPoints,
    const deformableBodyGeometricState& newStateDeformable,
    const bool exactDiameter
) const
{
    // TODO: find a software design approach to get rid of this dynamic cast!!
    const axisymmetricBodyGeometricState& newState 
        = dynamic_cast<const axisymmetricBodyGeometricState&>(newStateDeformable);

    // construct new point field centered at the origin
    pointField newPoints(undisplacedPoints - center());

    // scale points to new diameter, while conserving the volume
    tensor rotationToXAxis = 
                positiveDeterminantRotationTensor(axis(), vector(1,0,0));

    scalar diameterFactor;
    scalar axisFactor;
    if (exactDiameter)
    {
        diameterFactor = newState.diameter()/diameter();
        axisFactor = 1.0/pow(diameterFactor, 2);
    }
    else
    {
        axisFactor = newState.length()/length();
        diameterFactor = 1.0/sqrt(axisFactor);
    }
        
    tensor scaling(axisFactor,     0, 0,
                   0, diameterFactor, 0,
                   0, 0, diameterFactor);

    newPoints = inv(rotationToXAxis) & scaling & rotationToXAxis & newPoints;

    // displace points
    vector n1  = axis_/mag(axis_);
    vector n2  = newState.axis()/mag(newState.axis());
    tensor T   = positiveDeterminantRotationTensor(n1, n2);

    newPoints = (T & newPoints) + newState.center();

    return newPoints;
}


Foam::Istream&
Foam::axisymmetricBodyGeometricState::read(Istream& is)
{
    deformableBodyGeometricState::read(is);

    dictionary dict(is);
    axis_     = dict.lookup("axis");
    diameter_ = readScalar(dict.lookup("diameter"));
    length_   = readScalar(dict.lookup("length"));

    return is;
}


void
Foam::axisymmetricBodyGeometricState::write(Ostream& os) const
{
    deformableBodyGeometricState::write(os);
    os.writeKeyword("axis")     << axis_     << token::END_STATEMENT << nl;
    os.writeKeyword("diameter") << diameter_ << token::END_STATEMENT << nl;
    os.writeKeyword("length")   << length_   << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::axisymmetricBodyGeometricState::operator=(const axisymmetricBodyGeometricState& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::axisymmetricBodyGeometricState::operator="
                     "(const Foam::axisymmetricBodyGeometricState&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    deformableBodyGeometricState::operator=(rhs);

    axis_     = rhs.axis_;
    diameter_ = rhs.diameter_;
    length_   = rhs.length_;
}


// * * * * * * * * * * * * * * IOstream Operators * * * * * * * * * * * * * * //

Foam::Istream&
Foam::operator>>(Istream& is, axisymmetricBodyGeometricState& obj)
{
    obj.read(is);

    return is;
}


Foam::Ostream&
Foam::operator<<(Ostream& os, const axisymmetricBodyGeometricState& obj)
{
    obj.write(os);

    os.check("Ostream& operator<<(Ostream&, "
             "const axisymmetricBodyGeometricState&");

    return os;
}


// ************************************************************************* //
