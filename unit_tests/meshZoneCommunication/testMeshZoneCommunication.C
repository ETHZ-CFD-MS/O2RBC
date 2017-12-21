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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Output mesh information on terminal
    Pout<< "Number of points in mesh = " << mesh.nPoints() << endl;
    
    // Send mesh from master to processor 1
    
    runTime++;
    runTime.write();
    
    if (Pstream::master())
    {
        label destinationProc = 1;

        OPstream toNext
        (   
            Pstream::blocking,
            destinationProc
        );

        // Sending arrays
        const pointField& points = mesh.points();
        const faceList& faces = mesh.faces();
        const labelList& owner = mesh.faceOwner();
        const labelList& neighbour = mesh.faceNeighbour();

        toNext << points << faces << owner << neighbour;

        Pout<< "Points: " << points[0] << endl;
        Pout<< "Faces: " << faces[0] << endl;
        Pout<< "Owner: " << owner[0] << ", " << owner[1] << endl;
        Pout<< "Neighbour: " << neighbour[0] << ", " 
            << neighbour[1] << endl;
    }
    else if (Pstream::myProcNo() == 1)
    {
        label sourceProc = 0;

        IPstream fromPrevious
        (
            Pstream::blocking,
            sourceProc
        );

        pointField points;
        faceList   faces;
        labelList  owner;
        labelList  neighbour;

        fromPrevious >> points >> faces >> owner >> neighbour;

        Pout<< "points[0]:" << points[0] << endl;
        Pout<< "faces[0]: " << faces[0] << endl;
        Pout<< "owner: " << owner[0] << ", " << owner[1] << ", "
            << owner[2] << ", " << owner[3] << endl;
        Pout<< "neighbour: " << neighbour[0] << ", " << neighbour[1] << ", "
            << neighbour[2] << ", " << neighbour[3] << endl;

        // Create mesh
        Pout<< "Creating mesh" << endl;
        polyMesh newMesh
        (
            IOobject
            (
                "newMesh",
                mesh.time().constant(),
                mesh.time(),
                IOobject::NO_READ
            ),
            xferMove(points),
            xferMove(faces),
            xferMove(owner),
            xferMove(neighbour),
            false
        );
        newMesh.write();
        Pout<< "Created mesh" << endl;


    }

    // Output mesh information on terminal
    Pout<< "Number of points in mesh = " << mesh.nPoints() << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //




