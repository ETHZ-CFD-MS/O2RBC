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

Application
    testEllipticCylinder

Description
    Test the transformation of a straight cylinder to an elliptic cylinder
    given by vertex positions and associated radii.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "unitConversion.H"
#include "mathematicalConstants.H"
#include "tensor.H"
#include "transformField.H"

// return construct the shear transformation

// return bisection angle of a vertex
scalar bisectionAngle(const point& P, const point& Q, const point& R);

// return unit normal vector to vertex plane. If the points are colinear, an 
// arbitrary normal vector is returned.
vector normalVertexPlane(const point& P, const point& Q, const point& R);

// return unit bisection vector of a vertex
vector bisectionVector(const point& P, const point& Q, const point& R);

// return tensor that puts the elliptic cylinder along the x-axis, with the
// plane-normal vector at P along the z-axis
tensor transformToBase(const point& P, const point& Q, const vector nP);

// construct the tensor for a shear transformation
// planeNormal is the normal vector to the plane fixed by the shear
// transformation. imageNormal is the image of that vector under this
// transformation.
tensor shearTensor(const vector& planeNormal, const vector& imageNormal);

// compute the cylinder coordinates of a point (main axis = x-axis)
void cylinderCoords(const point& p, scalar& r, scalar& theta);

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // define vertices and their radii
    point vA(0,1,0);
    point vB(1,1,0);
    point vC(2,1,0);
    point vD(3,1,0);

    // scalar rA = 0.1;
    scalar rB = 0.1;
    scalar rC = 0.2;
    // scalar rD = 0.1;

    scalar lBC = mag(vC - vB);

    // get bisection angles
    scalar betaB = bisectionAngle(vA, vB, vC);
    scalar betaC = bisectionAngle(vB, vC, vD);

    Info << "Bisection angles: " << radToDeg(betaB) << ", " << radToDeg(betaC) << endl;

    // check if either angle is zero

    // get normal vector
    vector nB = normalVertexPlane(vA, vB, vC);
    vector nC = normalVertexPlane(vB, vC, vD);

    Info << "Normal vectors: " << nB << ", " << nC << endl;

    // get bisection vectors
    vector mB = bisectionVector(vA, vB, vC);
    vector mC = bisectionVector(vB, vC, vD);

    Info << "Bisection vectors: " << mB << ", " << mC << endl;

    // Now we have an orthogonal base at B and C:
    // (C-B, mB, nB) and (B-C, mC, nC)
    
    // Define a transformation that maps
    //    C-B ----> x-axis
    //     nB ----> e_z 
    tensor T    = transformToBase(vB, vC, nB);
    tensor Tinv = inv(T);

    vector TmB = transform(T, mB);
    vector TnB = transform(T, nB);
    vector TmC = transform(T, mC);
    vector TnC = transform(T, nC);

    Info << "Transformation tensor: " << T << endl;
    Info << "Image of vC - vB: " << transform(T, vC-vB) << endl;
    Info << "Image of mB:      " << TmB << endl;
    Info << "Image of nB:      " << TnB << endl;
    Info << "Image of mC:      " << TmC << endl;
    Info << "Image of nC:      " << TnC << endl;

    // get the length and the radius of the original mesh
    scalar L0 = mesh.bounds().span().x();
    scalar r0 = 0.5*mesh.bounds().span().y();
    point  pRight(L0, 0, 0);
    
    Info << "Original mesh size: L = " << L0 << ", r0 = " << r0 << endl;

    // find shear transformations for both sides of the cylinder
    tensor shearB = tensor::I;
    vector TmBproj(0.0, TmB.y(), TmB.z());
    if (mag(TmB - TmBproj) > SMALL)
    {
        shearB = shearTensor(TmBproj, TmB);
    }

    Info << "shearB = " << shearB << endl;
    
    tensor shearC = tensor::I;
    vector TmCproj(0.0, TmC.y(), TmC.z());
    if (mag(TmC - TmCproj) > SMALL)
    {
        shearC = shearTensor(TmCproj, TmC);
    }
    Info << "shearC = " << shearC << endl;

    // transform mesh
    pointField newPoints = mesh.points();

    forAll(mesh.points(), i)
    {
        const point p = mesh.points()[i];

        // construct scaled projections of p on cylinder sides
        vector xB = vector(0.0, rB/r0*p.y(), rB/r0*p.z());
        vector xC = vector(lBC, rC/r0*p.y(), rC/r0*p.z());

        // shear xB so that it is on the left side of the elliptic cylinder
        xB = transform(shearB, xB);

        // shear xC so that it is on the right side of the elliptic cylinder
        xC = transform(shearC, xC - pRight) + pRight;

        // compute convex combination of xB and xC
        scalar s = p.x()/L0;
        newPoints[i] = (1.0 - s)*xB + s*xC;

        // map the point to the graph edge position
        newPoints[i] = xB + (Tinv & newPoints[i]);
    }

    // write it
    runTime++;
    mesh.movePoints(newPoints);
    runTime.write();


    Info<< "End\n" << endl;

    return 0;
}

scalar bisectionAngle(const point& P, const point& Q, const point& R)
{
    vector QP = P - Q;
    vector QR = R - Q;

    // angle between both edges
    scalar alpha = Foam::acos((QP & QR)/(mag(QP)*mag(QR)));

    // bisection angle
    return 0.5*(constant::mathematical::pi-alpha);
}

vector normalVertexPlane(const point& P, const point& Q, const point& R)
{
    vector QP = P - Q;
    vector QR = R - Q;

    vector crossProduct = QP^QR;

    if (crossProduct == vector::zero)
    {
        // if QP and QR are colinear, just find any vector that is orthogonal.
        if (QP.x() != 0.0)
        {
            vector v(-QP.y()/QP.x(), 1.0, 0.0);
            return v/mag(v);
        }
        else if (QP.y() != 0.0)
        {
            vector v(1.0, -QP.x()/QP.y(), 0.0);
            return v/mag(v);
        }
        else if (QP.z() != 0.0)
        {
            vector v(0.0, 1.0, -QP.y()/QP.z());
            return v/mag(v);
        }
    }

    return crossProduct/mag(crossProduct);
}

vector bisectionVector(const point& P, const point& Q, const point& R)
{
    vector QP = P - Q;
    vector QR = R - Q;

    vector bisectV = QP/mag(QP) + QR/mag(QR);

    if (bisectV == vector::zero)
    {
        bisectV = QP ^ normalVertexPlane(P, Q, R);
    }

    return bisectV/mag(bisectV);
}

tensor transformToBase(const point& P, const point& Q, const vector nP)
{
    vector PQ = Q - P;
    // check that nP is orthogonal to PQ
    if (mag(nP & PQ) > SMALL)
    {
        FatalErrorIn
        (
            "transformToBase(const point&, const point&, const vector&)"
        ) << " The vector nP is not orthogonal to PQ." << nl
          << abort(FatalError);
    }

    // tensor that rotates PQ to x-axis
    vector ePQ = PQ/mag(PQ);
    tensor T1 = rotationTensor(ePQ, vector(1,0,0)); 

    // tensor that rotates T1(nP) to z-axis
    tensor T2 = rotationTensor(T1 & nP, vector(0,0,1));

    // return composition
    return T2 & T1;
}

void cylinderCoords(const point& p, scalar& r, scalar& theta)
{
    r = Foam::sqrt(p.y()*p.y() + p.z()*p.z());
    theta = Foam::atan2(p.z(), p.y());
}

tensor shearTensor(const vector& planeNormal, const vector& imageNormal)
{
    Info << "planeNormal = " << planeNormal << endl;
    Info << "imageNormal = " << imageNormal << endl;
    // check that imageNormal-planeNormal is orthogonal to planeNormal
    if (mag((imageNormal - planeNormal) & planeNormal) > SMALL)
    {
        FatalErrorIn
        (
            "tensor shearTensor(const vector&, const vector&)"
        ) << "The arguments do not define a shear transformation." << nl
          << "Check that imageNormal - planeNormal is orthogonal to planeNormal" << nl
          << abort(FatalError);
    }

    // construct transformation so that shear fixes (x,z)-plane and is along
    // x-direction:
    // the tensor T1 maps the plane-normal vector to the y-axis
    tensor T1 = rotationTensor(planeNormal/mag(planeNormal), vector(0,1,0));
    vector nT1diff = transform(T1, imageNormal-planeNormal);
    nT1diff = nT1diff/mag(nT1diff);
    // the tensor T2 then rotates the shearing vector (i.e., the difference
    // between imageNormal and planeNormal) so it is along the x-axis.
    tensor T2 = rotationTensor(nT1diff, vector(1,0,0));
    // Info << T1 << nl << T2 << endl;

    tensor T_rot = transform(T2, T1);

    // check that T2(T1(im)) has same y-coordinate as T2(T1(pl))
    vector T2T1pl = transform(T_rot, planeNormal);
    vector T2T1im = transform(T_rot, imageNormal);
    if (mag(T2T1pl.y() - T2T1im.y()) > SMALL)
    {
        Info << "ERROR, check code..." << endl;
    }

    scalar lambda = T2T1im.x() / T2T1im.y();
    Info << "lambda = " << lambda << endl;

    // Construct the shear tensor (uses row-major ordering)
    tensor T_shear(1, lambda, 0, 0, 1, 0, 0, 0, 1);

    tensor result = inv(T_rot) & T_shear & T_rot;

    if (mag(transform(result, planeNormal) - imageNormal) > SMALL)
    {
        Info << "Wrong shear transformation..." << endl;
    }

    return result;
}
