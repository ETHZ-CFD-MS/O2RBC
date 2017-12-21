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

InNamespace
    Foam

Description
    Various functions to operate on points, vectors and tensors.

SourceFiles
    geometricOps.C

\*---------------------------------------------------------------------------*/

#ifndef geometricOps_H
#define geometricOps_H

#include "point.H"
#include "vector.H"
#include "tensor.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Return an arbitrarily constructed normal vector to the given vector.
vector normalVector(const vector& v);

//- Return the normal vector to both given vectors.
//  If the vectors are colinear, it returns an arbitrarily constructed vector
//  that is normal to them.
vector normalVector(const vector& v1, const vector& v2);

//- Return the angle between two vector
scalar angle(const vector& v1, const vector& v2);

//- Return the bisection vector formed by two vectors.
//  Returns a unit vector that forms an angle <= pi/2 with v1 and v2.
vector bisectionVector(const vector& v1, const vector& v2);

//- Return the major bisection vector formed by two vectors.
//  Returns a unit vector that forms an angle >= pi/4 with $\pm$ v1 and $\pm$ v2.
//  In this function, v1 and v2 are interpreted as lines rather than vector.
//  The returned vector forms the largest bisection angle between itself and
//  both $\pm$ v1 and $\pm$ v2.
vector majorBisectionVector(const vector& v1, const vector& v2);

//- Return the bisection angle formed by two vectors.
//  The returned angle is between 0 and pi/2.
scalar bisectionAngle(const vector& v1, const vector& v2);

//- Return the major bisection angle formed by two vectors.
//  The returned angle is between pi/4 and 3*pi/4.
scalar majorBisectionAngle(const vector& v1, const vector& v2);

//- Construct a rotation tensor that maps a pair of orthogonal vectors to
//  another pair of orthogonal vectors.
tensor rotationTensor(const Pair<vector>&, const Pair<vector>&);

//- Construct a rotation tensor that maps a unit vector to another, while avoiding
//  to create a tensor with determinant -1 (unlike the version in transform.H).
tensor positiveDeterminantRotationTensor(const vector& n1, const vector& n2);

//- Construct the tensor for a shear transformation.
//  The argument planeNormal is the normal vector to the plane fixed by the 
//  shear transformation. imageNormal is the image of that vector under this
//  transformation.
tensor shearTensor(const vector& planeNormal, const vector& imageNormal);

//- Return whether two vectors are colinear with a given angle tolerance
bool almostColinear(const vector& v1, const vector& v2, const scalar tol = 1e-6);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
