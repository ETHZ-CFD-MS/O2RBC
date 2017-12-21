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
    initializeHbPOGraph

Description
    Initialize the fields PO2 and Hb, and write them to the disk.
    
    PO2 is first read from the disk. 

    Hb is computed to be in equilibrium with PO2_RBC (with is interpolated from
    PO2).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "circularTubeGraph.H"

#include "RBCCollection.H"
#include "RBCInterpolation.H"
#include "RBCPathMover.H"
#include "dissociationCurve.H"

#include "unistd.h"
#include "stdio.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // For parallel debugging
    // {
        // int i = 0;
        // char hostname[256];
        // gethostname(hostname, sizeof(hostname));
        // Pout << "PID " << getpid() << " on " << hostname << " ready for attach " << endl;
        // while (0 == i)
            // Foam::sleep(5);
    // }

    //- Create
    //  * PO2, PO2Adv
    //  * PO2Plasma, PO2_RBC_inlet
    #include "createFields.H"

    // create dissociation curve
    autoPtr<dissociationCurve> DCPtr
    (
        dissociationCurve::New(transportProperties.subDict("dissociationCurve"))
    );

    // create RBCs
    RBCCollection RBCs(mesh, DCPtr());

    // create graph
    circularTubeGraph graph
    (
        Foam::IOobject
        (
            "graphDict",
            runTime.caseConstant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // create the RBC mover object
    RBCPathMover mover(mesh, RBCs);

    // set the initial position of the RBCs.
    mover.setInitialPositions();

    // create the interpolation object
    RBCInterpolation interpolator(mesh, RBCs);

    // interpolate PO2_RBC to PO2
    interpolator.interpolateRBCToEuler("PO2_RBC", PO2);

    // copy PO2 to PO2Adv
    PO2Adv = PO2;

    // write the initial positions of the meshes and the RBCs
    RBCs.RBCMesh().write();
    RBCs.write();

    // write fields at the initial time
    PO2.write();
    PO2Adv.write();
    RBCs.Hb().write();
    RBCs.PO2().write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
