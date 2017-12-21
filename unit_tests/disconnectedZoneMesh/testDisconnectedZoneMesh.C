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

#include "disconnectedZoneMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    disconnectedZoneMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    volScalarField c
    (
        IOobject
        (
            "c",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh 
    );

    Info<< "Cell zone names" << mesh.cellZones().names() << endl;
    Info<< "Point zone names" << mesh.pointZones().names() << endl;
    Info<< "Face zone names" << mesh.faceZones().names() << endl;

    Info<< "Bounds of cell zone 0: " << mesh.zoneBounds(0) << endl;
    Info<< "Bounds of cell zone 1: " << mesh.zoneBounds(1) << endl;

    Info<< "Field on zone 0: " << mesh.getFieldOnZone(c, 0) << endl;
    Info<< "Field on zone 1: " << mesh.getFieldOnZone(c, 1) << endl;

    Info<< "Boundary faces on zone 0: " << mesh.zoneBoundaryFaces(0) << endl;
    Info<< "Boundary faces on zone 1: " << mesh.zoneBoundaryFaces(1) << endl;

    mesh.setFieldOnZone(c, 2., 0);

    Info<< "Field after modification on zone 0: \n" << c << endl;

    mesh.setFieldOnZone(c, 3., 1);

    Info<< "Field after modification on zone 1: \n" << c << endl;

    c[0] = 123.;
    mesh.setBoundaryFieldToPatchInternalFieldOnZone(c, 0);
    Info<< "Field after boundary field modification: \n" << c << endl;

    Info<< "Mesh has non empty zone 0: " << mesh.hasNonEmptyZone(0) << endl;
    Info<< "Mesh has non empty zone 2: " << mesh.hasNonEmptyZone(2) << endl;


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

