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

#include "cartesianProcMeshInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    cartesianProcMeshInfo procMeshInfo(mesh);

    List<boundBox> procBBoxes = procMeshInfo.procBBoxes();
    forAll(procBBoxes, i)
    {
        Pout<< "Bounding box of processor mesh " << i
            << ": " << procBBoxes[i] << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

