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
    testDivergenceFreeProj

Description
    Test for projection of a non-divergence-free velocity field.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"

    solve
    (
        fvm::laplacian(f) == fvc::div(phi)
    );

    runTime++;

    surfaceScalarField Gf = fvc::snGrad(f) * mesh.magSf();
    psi = phi - Gf;


    phi.write();
    f.write();
    psi.write();

    Info<< "Maximum of divergence of psi: " << endl;
    Info<< max(mag(fvc::div(psi))) << endl;

    // V = linearInterpolate(psi);

    return 0;
}
