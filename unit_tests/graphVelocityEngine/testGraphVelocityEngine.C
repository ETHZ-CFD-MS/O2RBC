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
    testGraphVelocityEngine

Description
    Unit test for generating a divergence-free velocity field from velocities
    defined on a graph.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "circularTubeGraph.H"
#include "graphVelocityEngine.H"
#include "vascularGraphRegions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Create graph
    circularTubeGraph graph
    (
        Foam::IOobject
        (
            "graphDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read the regions for the vasculature
    vascularGraphRegions vesselRegions(mesh, graph);

    // Create the velocity engine
    graphVelocityEngine velocityEngine(mesh, "lumen", graph, vesselRegions);

    while(runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        const surfaceScalarField& conservativePhi = velocityEngine.conservativePhi();
        conservativePhi.write();

        runTime.write();

        Info<< endl;
    }


    Info<< "End\n" << endl;

    return 0;
}
