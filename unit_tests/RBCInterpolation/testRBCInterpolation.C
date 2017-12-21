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
    testRBCInterpolation

Description
    Unit test for interpolation between moving meshes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOmanip.H"

#include "cylinderGraph.H"
#include "RBC.H"
#include "regionProperties.H"
#include "dissociationCurve.H"
#include "RBCCollection.H"
#include "RBCPositionTester.H"
#include "RBCInterpolation.H"

#include "indexedOctree.H"
#include "treeDataCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void testVolume(const fvMesh& mesh, const RBCCollection& RBCs, const volScalarField& in_RBC_euler);


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    //- Create
    //  * PO2_crit, RBC_velocity
    //  * O2_consumption_rate
    //  * PO2_RBC_inlet, PO2_plasma, U_plasma
    //  * dissociation_rate, P_50, hill_exponent
    //  * kappa_Hb, alpha_RBC, vol_mol_O2, N_Hb
    //  * U, kappa_O2, in_tissue
    //  * PO2, PO2Adv, PO2Phi, Hb_euler, 
    //  * in_RBC_euler, in_RBC_mask
    //  * Mp, Mc, Rp, Rc
    //  * phi
    #include "createFields.H"

    regionProperties rp(runTime);

    // create dissociation curve
    dissociationCurve HbPO2_DC(P_50, hill_exponent);

    // create RBCs
    RBCCollection RBCs(runTime, rp.regionNames().size());

    forAll(rp.regionNames(), i)
    {
        Info << "Creating RBC " << rp.regionNames()[i] << nl << endl;
        RBCs.set
        (
            i,
            new RBC
            (
               rp.regionNames()[i],
               runTime,
               HbPO2_DC
            )
        );
    }

    RBCs[0].setCenter(point(12e-3,0,0));
    RBCs[0].RBCMesh().write();

    Info<< "Mesh bbox = " << mesh.bounds() << endl;
    Info<< "RBC bbox  = " << RBCs[0].RBCMesh().bounds() << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    // create the position helper
    RBCPositionTester positionTester(mesh, RBCs);

    // create the interpolation object
    RBCInterpolation interpolator(mesh, RBCs, positionTester);

    Info<< "Target to source weights: " << endl;
    Info<< interpolator.eulerToRBCInterpolator(0).tgtToSrcCellWght();
    Info<< "Source to target weights: " << endl;
    Info<< interpolator.eulerToRBCInterpolator(0).srcToTgtCellWght();

    // initialize fields from the "base" fields PO2, Hb and kappa_O2_no_RBC.
    #include "initializeFields.H"

    // Compute in_RBC_euler
    in_RBC_euler = 0;
    interpolator.interpolateRBCToEuler("in_RBC", in_RBC_euler);

    testVolume(mesh, RBCs, in_RBC_euler);

    // Move RBCs
    RBCs[0].setCenter(point(15e-3,0,0));

    // Update interpolator
    interpolator.update();

    Info<< "Target to source weights: " << endl;
    Info<< interpolator.eulerToRBCInterpolator(0).tgtToSrcCellWght();
    Info<< "Source to target weights: " << endl;
    Info<< interpolator.eulerToRBCInterpolator(0).srcToTgtCellWght();

    // Recompute in_RBC_euler
    in_RBC_euler = 0;
    interpolator.interpolateRBCToEuler("in_RBC", in_RBC_euler);

    testVolume(mesh, RBCs, in_RBC_euler);

    runTime++;
    runTime.write();

    // test findInside
    pointField points = RBCs[0].RBCMesh().points();
    forAll(points, i)
    {
        point pt = points[i];
        label cellIdx = mesh.cellTree().findInside(0.9*pt);
        Info << "Point " << pt << " found in cell with idx " << cellIdx << endl;
    }
    Info << "Mesh bounds = " << mesh.bounds() << endl;

    Info<< "End\n" << endl;

    return 0;
}

void testVolume(const fvMesh& mesh, const RBCCollection& RBCs, const volScalarField& in_RBC_euler)
{
    scalar RBCVol = RBCs[0].volume();
    scalar sumInRBCEuler = sum(in_RBC_euler*mesh.V()).value();
    Info<< "RBC volume: " << RBCVol << endl;
    Info<< "Weighted sum of in_RBC_euler volume: " << sumInRBCEuler<< endl;

    if (RBCVol != sumInRBCEuler)
    {
        Info<< "TEST FAILURE: Test of volume equivalence failed" << endl;
    }
}


// ************************************************************************* //


