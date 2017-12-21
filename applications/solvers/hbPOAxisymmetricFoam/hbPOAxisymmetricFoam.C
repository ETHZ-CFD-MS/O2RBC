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
    hbPOAxisymmetricFoam

Description
    Transient solver for advection-diffusion-reaction of PO2 with hemoglobin
    in discrete red blood cells (RBC) that flow through an axisymmetric domain.
    The RBC motion is dealt by an instance of RBCVelocityMover.

    The frame of reference of the computation is Eulerian.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOmanip.H"
#include "expressionSource.H"
#include "processorFvPatch.H"

#include "RBCCollection.H"
#include "RBCInterpolation.H"
#include "RBCVelocityMover.H"
#include "dissociationCurve.H"

#include "implicitSourceControl.H"

#include <boost/math/special_functions/fpclassify.hpp>
// For parallel debugging
#include "unistd.h"
#include "stdio.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void computePhi
(
    surfaceScalarField& phi,
    const scalar RBCVelocity, 
    const fvMesh& mesh,
    const word& plasmaZoneName
);

bool catchSmallOrNanValues(const volScalarField&);

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
    
    //- Create fields
    #include "createFields.H"

    const word plasmaZoneName("plasma");

    // Create dissociation curve
    autoPtr<dissociationCurve> DCPtr
    (
        dissociationCurve::New(transportProperties.subDict("dissociationCurve"))
    );
    const dissociationCurve& DC = DCPtr();

    // Create RBCs
    RBCCollection RBCs(mesh, DC);

    // reaction factor in the O2 diffusion-reaction equation
    dimensionedScalar O2ReactionFactor = NHb*VMolO2*dissociationRate; // [mlO2 m^-3 s^-1]

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    implicitSourceControl implicitSource(mesh);

    // create the interpolation object
    RBCInterpolation interpolator(mesh, RBCs);

    // create the RBC mover object
    autoPtr<RBCVelocityMover> moverPtr(RBCVelocityMover::New(mesh, RBCs));

    // initialize dependent fields from base fields, and write them to disk
    #include "initializeFields.H"

    #include "computeXiRBC.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // Move RBCs
        moverPtr->moveAll();

        // Compute the advective flux
        computePhi(phi, moverPtr->getRBCVelocity(), mesh, plasmaZoneName);

        // Update interpolators
        interpolator.update();

        // Recompute fields that depend on the RBC position
        oldInRBCEuler.field() = inRBCEuler.field();
        inRBCEuler = 0;
        inRBCMask  = 0;
        interpolator.interpolateRBCToEuler("in_RBC", inRBCEuler);

        forAll(mesh.C(), cI)
        {
            if (inRBCEuler[cI] > 0)
            {
                inRBCMask[cI] = 1;
            }
        }

        volScalarField alphaPrevious("alphaPrevious", alpha);
        alpha   = inRBCEuler*alphaRBC   + (1 - inRBCEuler)*alphaNoRBC;
        kappaO2 = inRBCEuler*kappaO2RBC + (1 - inRBCEuler)*kappaO2NoRBC;

        #include "PO2AdvEqn.H"

        while (implicitSource.loop())
        {
            // interpolate fields from the previous iteration or time step
            HbEuler = 0.0;
            interpolator.interpolateRBCToEuler("Hb", HbEuler);
            #include "PO2Eqn.H"
            interpolator.interpolateEulerToRBC(PO2, "PO2");
            #include "HbEqn.H"
        }
        
        // update oxygen concentration
        C = alpha*PO2;

        O2ConsumptionRateOutput = O2ConsumptionRate();
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //


void computePhi
(
    surfaceScalarField& phi,
    const scalar RBCVelocity, 
    const fvMesh& mesh,
    const word& plasmaZoneName
)
{
    const cellZoneMesh& czMesh = mesh.cellZones();
    const cellZone& plasmaZone = czMesh[plasmaZoneName];
    volVectorField U
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("U", dimVelocity, vector::zero)
    );
    // set values of U in the domain interior
    forAll(plasmaZone, cI)
    {
        U[cI] = vector(RBCVelocity, 0.0, 0.0);
    }
    // set values of U on patches
    forAll(mesh.boundaryMesh(), patchI)
    {
        U.boundaryFieldRef()[patchI] = U.boundaryField()[patchI].patchInternalField();
    }
    phi = linearInterpolate(U) & mesh.Sf();
}


bool catchSmallOrNanValues(const volScalarField& f)
{
    bool caughtValue = false;

    if (isnan(sum(f).value()) || min(f).value() < 1e-8)
    {
        Pout<< "Unexpected value of field " << f.name() << endl;
        Pout<< "min = " << min(f).value() << ", "
            << "max = " << max(f).value() << endl;

        caughtValue = true;

        label counter = 0;
        forAll(f, i)
        {
            if ((f[i] != f[i] || f[i] < 1.0e-1) && counter <= 10)
            {
                Pout<< "For cell " << i << ", " << f.name() 
                    << " = " << f[i]
                    << "(faces = " << f.mesh().cells()[i] << ")" << endl;
                counter++;
            }
        }

        forAll(f.mesh().boundary(), patchI)
        {
            label counter = 0;
            forAll(f.boundaryField()[patchI], i)
            {
                scalar val = f.boundaryField()[patchI][i];
                if ((val != val || val < 1.0e-1) && counter <= 10)
                {
                    Pout<< "For boundary face " << i << ", " << f.name()
                        << "= " << val << endl;
                    counter++;
                }
            }
        }
    }

    return caughtValue;
}


