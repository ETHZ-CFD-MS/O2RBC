/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    adjustRBCMeshVolume

Description
    Adjust the volume of a RBC mesh by rescaling it.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "deformableBodyGeometricState.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Adjust the volume of a RBC mesh"
    );
    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    #include "createNamedMesh.H"
    const word oldInstance = mesh.pointsInstance();

    scalar originalVolume = sum(mesh.V()).value();

    IOdictionary dict
    (
        IOobject
        (
            "RBCBaseGeometricState",
            runTime.caseConstant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false // do not register in database
        )
    );

    autoPtr<deformableBodyGeometricState> pGeometricState 
        = deformableBodyGeometricState::New(dict, mesh);
    vector center = pGeometricState->center();
    scalar prescribedVolume = readScalar(dict.lookup("volume"));

    pointField pf = mesh.points();
    pf -= center;
    pf *= std::pow(prescribedVolume/originalVolume, 1./3.);
    pf += center;
    mesh.movePoints(pf);

    mesh.setInstance(oldInstance);
    mesh.write();

    Info<< "RBC mesh volume adjusted from " << originalVolume
        << " to " << sum(mesh.V()).value() << "." << nl
        << "Bounding box: " << mesh.bounds() << nl
        << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

