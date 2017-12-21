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
    hbPOGraphFoam

Description
    Transient solver for the advection-diffusion-reaction of PO2 and hemoglobin 
    saturation in discrete moving red blood cells (RBC). Uses provided paths 
    for RBCs through a graph structure.

    The frame of reference of the computation is Eulerian.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOmanip.H"
#include "expressionSource.H"

#include "circularTubeGraph.H"
#include "vascularGraphRegions.H"
#include "regionDependentField.H"
#include "graphVelocityEngine.H"

#include "RBCCollection.H"
#include "RBCInterpolation.H"
#include "RBCPathMover.H"
#include "dissociationCurve.H"

#include "implicitSourceControl.H"

#include <boost/math/special_functions/fpclassify.hpp>
// For parallel debugging
#include "unistd.h"
#include "stdio.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
    
    // Create graph
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

    // Read the regions for the vasculature
    vascularGraphRegions vesselRegions(mesh, graph);

    // Create the velocity engine
    graphVelocityEngine velocityEngine(mesh, "lumen", graph, vesselRegions);

    //- Create fields
    //  * PO2Crit, O2ConsumptionRate
    //  * dissociationRate
    //  * kappaHb, VMolO2, NHb
    //  * PO2, PO2Adv, PO2Phi, C, Hb_euler, 
    //  * in_RBC_euler, in_RBC_mask
    #include "createFields.H"

    // Create dissociation curve
    autoPtr<dissociationCurve> DCPtr
    (
        dissociationCurve::New(transportProperties.subDict("dissociationCurve"))
    );
    const dissociationCurve& DC = DCPtr();

    // Create RBCs
    RBCCollection RBCs(mesh, DC);

    // reaction factor in the O2 diffusion-reaction equation
    dimensionedScalar O2_reaction_factor = NHb*VMolO2*dissociationRate; // [mlO2 m^-3 s^-1]

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    implicitSourceControl implicitSource(mesh);

    // create the interpolation object
    RBCInterpolation interpolator(mesh, RBCs);

    // create the RBC mover object
    RBCPathMover mover(mesh, RBCs);

    // initialize fields from the "base" fields PO2, Hb and kappa_O2_no_RBC.
    #include "initializeFields.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // Move RBCs
        mover.moveAll();

        // Compute the advective flux
        surfaceScalarField phi = velocityEngine.conservativePhi();

        // Compute the upwind cells
        // #include "computeUpwindCells.H"

        // Update interpolators
        interpolator.update();

        // Recompute in_RBC_euler, in_RBC_mask and kappa_O2
        in_RBC_euler = 0;
        in_RBC_mask  = 0;
        interpolator.interpolateRBCToEuler("in_RBC", in_RBC_euler);

        forAll(mesh.C(), cI)
        {
            if (in_RBC_euler[cI] > 0)
            {
                in_RBC_mask[cI] = 1;
            }
        }

        volScalarField alphaPrevious = alpha.field();
        alpha.update(in_RBC_euler);
        kappa_O2.update(in_RBC_euler);

        #include "PO2AdvEqn.H"

        while (implicitSource.loop())
        {
            // interpolate fields from the previous iteration or time step
            Hb_euler = 0.0;
            interpolator.interpolateRBCToEuler("Hb", Hb_euler);
            interpolator.interpolateEulerToRBC(PO2, "PO2");

            #include "PO2Eqn.H"
            #include "HbEqn.H"
        }
        
        // update oxygen concentration
        C = alpha.field()*PO2;

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


