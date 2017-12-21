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
    testMeshZoneCommunication

Description
    Unit test for sending mesh zones across processors.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOmanip.H"

#include "regionAddRemoveFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void displayMeshInfo(const fvMesh& mesh);

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    regionAddRemoveFvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    polyMesh RBCRegion
    (
        IOobject
        (
            "RBC2",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    volScalarField f
    (
        IOobject
        (
            "f",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    );

    forAll(mesh.cells(), i)
    {
        f[i] = i;
    }
    f.write();


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    displayMeshInfo(mesh);

    autoPtr<mapPolyMesh> meshMap;

    runTime++;
    mesh.removeMesh("RBC1", "RBC1Points", "RBC1Faces");
    mesh.update();
    displayMeshInfo(mesh);
    mesh.write();

    runTime++;
    mesh.addMesh(RBCRegion, true);
    mesh.merge();
    mesh.update();
    displayMeshInfo(mesh);
    mesh.write();

    // runTime++;
    // mesh.removeMesh("RBC", "RBCPoints", "RBCFaces");
    // mesh.update();
    // displayMeshInfo(mesh);
    // mesh.write();

    runTime++;
    mesh.removeMesh("RBC", "RBCPoints", "RBCFaces");
    mesh.update();
    displayMeshInfo(mesh);
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //


void displayMeshInfo(const fvMesh& mesh)
{
    Info<< "Mesh information at time " << mesh.time().timeName() << nl
        << "   Number of cells: " << mesh.nCells() << nl
        << "   Number of points: " << mesh.nPoints() << nl
        << "   Patch names: " << mesh.boundaryMesh().names() << nl
        << "   Cell zones: " << mesh.cellZones().names() << nl
        << "   Point zones: " << mesh.pointZones().names() << nl
        << "   Face zones: " << mesh.faceZones().names() << nl
        << endl;
}


