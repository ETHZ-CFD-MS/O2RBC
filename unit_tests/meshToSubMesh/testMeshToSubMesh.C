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
    testMeshToSubMesh

Description
    Unit test for copying fields between a mesh and a submesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "meshToSubMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    fvMesh subMesh
    (
        IOobject
        (
            "subMesh",
            mesh.time().timeName(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    volScalarField f
    (
        IOobject
        (
            "f",
            mesh.time().timeName(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField fSub
    (
        IOobject
        (
            "fSub",
            subMesh.time().timeName(),
            subMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        subMesh,
        f.dimensions()
    );

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    surfaceScalarField phiSub
    (
        IOobject
        (
            "phiSub",
            subMesh.time().timeName(),
            subMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        subMesh,
        phi.dimensions()
    );
    phiSub = 0.0;

    meshToSubMesh copyObj(mesh, subMesh);

    // Test the copy to subMesh
    copyObj.copyMeshToSubMesh(f, fSub);
    copyObj.copyMeshToSubMesh(phi, phiSub);

    fSub.write();
    phiSub.write();

    // Test the copy to the mesh
    runTime++;
    fSub   = 2;
    phiSub = 2;
    copyObj.copySubMeshToMesh(f, fSub);
    copyObj.copySubMeshToMesh(phi, phiSub);

    runTime.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

