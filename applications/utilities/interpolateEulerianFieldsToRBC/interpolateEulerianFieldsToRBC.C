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
    interpolateEulerianFieldsToRBC

Description
    Interpolate the PO2 field from the Eulerian to the RBC mesh. Then, set
    the hemoglobin saturation in the RBC to be in equilibrium with PO2.

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

    //- Create
    //  * transportProperties, PO2
    #include "createFields.H"

    autoPtr<dissociationCurve> DCPtr
    (
        dissociationCurve::New(transportProperties.subDict("dissociationCurve"))
    );

    RBCCollection RBCs(mesh, DCPtr());

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    RBCInterpolation interpolator(mesh, RBCs);

    interpolator.interpolateEulerToRBC(PO2, "PO2");
    RBCs.setEquilibriumHb();

    RBCs.Hb().write();
    RBCs.PO2().write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
