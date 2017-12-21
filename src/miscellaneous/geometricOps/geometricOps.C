/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "geometricOps.H"

#include "mathematicalConstants.H"
#include "transform.H"


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

Foam::vector Foam::normalVector(const vector& v)
{
    vector n = vector::zero;

    if (v.x() != 0.0)
    {
        n = vector(-v.y()/v.x(), 1.0, 0.0);
    }
    else if (v.y() != 0.0)
    {
        n = vector(1.0, -v.x()/v.y(), 0.0);
    }
    else if (v.z() != 0.0)
    {
        n = vector(0.0, 1.0, -v.y()/v.z());
    }
    else
    {
        FatalErrorIn
        (
            "Foam::normalVector(const vector& v) const"
        ) << "The normal vector to the zero vector is not well defined." << nl
          << abort(FatalError);
    }

    return n/mag(n);
}


Foam::vector Foam::normalVector
(
    const vector& v1,
    const vector& v2
)
{
    vector normal = vector::zero;
    vector crossProduct = v1^v2;

    if (almostColinear(v1, v2))
    {
        normal = normalVector(v1);
    }
    else
    {
        normal = crossProduct;
    }

    return normal/mag(normal);
}


Foam::scalar Foam::angle
(
    const vector& v1,
    const vector& v2
)
{
    scalar cosAlpha = (v1 & v2)/(mag(v1)*mag(v2));

    // correct cosAlpha in case numerical errors cause values slightly larger
    // than 1 or smaller than -1.
    if (cosAlpha > 1 && cosAlpha < 1 + SMALL)
    {
        cosAlpha = 1;
    }
    else if (cosAlpha < -1 && cosAlpha > -1 - SMALL)
    {
        cosAlpha = -1;
    }

    // the angle alpha takes a value between 0 and pi.
    scalar alpha = Foam::acos(cosAlpha);

    return alpha;
}


Foam::vector Foam::bisectionVector
(
    const vector& v1,
    const vector& v2
)
{
    vector bisect = v1/mag(v1) + v2/mag(v2);

    if (almostColinear(v1, v2))
    {
        bisect = v1 ^ normalVector(v1);
    }

    return bisect/mag(bisect);
}


Foam::vector Foam::majorBisectionVector
(
    const vector& v1,
    const vector& v2
)
{
    if (angle(v1,v2) >= 0.5*constant::mathematical::pi)
    {
        return bisectionVector(v1,v2);
    }
    else
    {
        return bisectionVector(v1,-v2);
    }
}


Foam::scalar Foam::bisectionAngle
(
    const vector& v1,
    const vector& v2
)
{
    return 0.5*angle(v1,v2);
}


Foam::scalar Foam::majorBisectionAngle
(
    const vector& v1,
    const vector& v2
)
{
    if (angle(v1,v2) >= 0.5*constant::mathematical::pi)
    {
        return 0.5*angle(v1,v2);
    }
    else
    {
        return 0.5*(constant::mathematical::pi - angle(v1,v2));
    }
}


Foam::tensor
Foam::rotationTensor(const Pair<vector>& vv, const Pair<vector>& ww)
{
    // check orthogonality of argument vectors
    if (abs(vv.first() & vv.second()) > SMALL ||
        abs(ww.first() & ww.second()) > SMALL)
    {
        FatalErrorIn("Foam::rotationTensor(const Pair<vector>&, const Pair<vector>&)")
            << "The pair of vectors passed as arguments are not orthogonal." << nl
            << "vv = " << vv << ", ww = " << ww << nl
            << abort(FatalError);
    }

    const vector v1 = vv.first()/mag(vv.first());
    const vector v2 = vv.second()/mag(vv.second());
    const vector v3 = v1^v2;

    const vector w1 = ww.first()/mag(ww.first());
    const vector w2 = ww.second()/mag(ww.second());
    const vector w3 = w1^w2;

    tensor Tv = (tensor(v1, v2, v3)).T();
    tensor Tw = (tensor(w1, w2, w3)).T();

    tensor result = Tw & inv(Tv);

    return result;
}


Foam::tensor 
Foam::positiveDeterminantRotationTensor(const vector& n1, const vector& n2)
{
    const scalar s = n1 & n2;

    if (s < -1 + SMALL)
    {
        const vector n3 = normalVector(n1);
        return rotationTensor(Pair<vector>(n1, n3), Pair<vector>(n2, n3));
    }
    else
    {
        return rotationTensor(n1, n2);
    }
}

Foam::tensor 
Foam::shearTensor
(
    const vector& planeNormal,
    const vector& imageNormal
)
{
    // if imageNormal is equal to planeNormal, return identity tensor
    if (mag(imageNormal - planeNormal) < SMALL)
    {
        return tensor::I;
    }

    // check that imageNormal-planeNormal is orthogonal to planeNormal
    if (mag((imageNormal - planeNormal) & planeNormal) > SMALL)
    {
        FatalErrorIn
        (
            "Foam::shearTensor(const vector&, const vector&)"
        ) << "The arguments do not define a shear transformation." << nl
          << "Check that (imageNormal - planeNormal) is orthogonal to "
          << "planeNormal" << nl
          << abort(FatalError);
    }

    // construct tensor to rotate (shearV, planeNormal) to (e_x, e_y)
    vector unitPlaneNormal = planeNormal/mag(planeNormal);
    vector unitShearVector = (imageNormal - planeNormal)
                            / mag(imageNormal - planeNormal);
    tensor T = rotationTensor(Pair<vector>(unitShearVector, unitPlaneNormal),
                              Pair<vector>(vector(1,0,0), vector(0,1,0)));

    vector Tpl = T & planeNormal;
    vector Tim = T & imageNormal;

    // check that T(imageNormal)) has same y-coordinate as T(planeNormal))
    if (mag(Tpl.y() - Tim.y()) > SMALL)
    {
        FatalErrorIn
        (
            "Foam::shearTensor(const vector&, const vector&)"
        ) << "The image of planeNormal and imageNormal do not have the same" << nl
          << "y-coordinate. There is a bug in the code." << nl
          << abort(FatalError);
    }

    scalar lambda = Tim.x()/Tim.y();

    // Construct the shear tensor (uses row-major ordering)
    tensor T_shear(1, lambda, 0, 0, 1, 0, 0, 0, 1);

    tensor result = inv(T) & T_shear & T;

    // sanity check
    if (mag((result & planeNormal) - imageNormal) > SMALL)
    {
        FatalErrorIn
        (
            "Foam::shearTensor(const vector&, const vector&)"
        ) << "The constructed tensor does not map planeNormal onto imageNormal." << nl
          << "y-coordinate. There is a bug in the code." << nl
          << abort(FatalError);
    }

    return result;
}


bool Foam::almostColinear
(
    const vector& v1,
    const vector& v2,
    const scalar tol
)
{
    scalar ang = angle(v1, v2);
    if (mag(ang) < tol || mag(constant::mathematical::pi - ang) < tol)
    {
        return true;
    }
    else
    {
        return false;
    }
}
