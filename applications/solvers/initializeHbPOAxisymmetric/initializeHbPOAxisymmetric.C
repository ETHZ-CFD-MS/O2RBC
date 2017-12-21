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
    initializeHbPOEuler

Description
    Initialize the fields PO2 and Hb, and write them to the disk.
    
    PO2 is first read from the disk. In the Lagrangian meshes, PO2_lag is read and
    is interpolated to PO2 and PO2_RBC.

    Hb is computed to be in equilibrium with PO2_RBC (with is interpolated from
    PO2_lag).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "RBCCollection.H"
#include "RBCInterpolation.H"
#include "RBCVelocityMover.H"
#include "dissociationCurve.H"

#include "analyticalSolutionPO2Cone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"

    // create dissociation curve
    autoPtr<dissociationCurve> DCPtr
    (
        dissociationCurve::New(transportProperties.subDict("dissociationCurve"))
    );

    // create RBCs
    RBCCollection RBCs(mesh, DCPtr());

    // create analytical solver
    analyticalSolutionPO2Cone analyticalSolver(transportProperties, geometricProperties);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    RBCVelocityMover mover(mesh, RBCs);

    mover.setInitialPositions();

    RBCInterpolation interpolator(mesh, RBCs);

    analyticalSolver.computePO2(PO2);
    analyticalSolver.computeHbAndRBCPO2(RBCs);

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
