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
    setSurfaceDistanceField

Description
    Sets value of a field based on the distance to a STL surface.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "triSurfaceDistanceFunction.H"

int main(int argc, char *argv[])
{
    timeSelector::addOptions(true, false);
    argList::addOption
    (
        "dict",
        "file",
        "specify an alternative dictionary for the setSurfaceDistanceField dictionary"
    );

    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    const word dictName("setSurfaceDistanceFieldDict");

    fileName dictPath = dictName;
    if (args.optionFound("dict"))
    {
        dictPath = args["dict"];
        if (isDir(dictPath))
        {
            dictPath = dictPath / dictName;
        }
    }

    Info<< "Reading " << dictName << "\n" << endl;

    IOdictionary setSurfaceDistanceFieldDict
    (
        (
            args.optionFound("dict")
          ? IOobject
            (
                dictPath,
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
          : IOobject
            (
                dictName,
                runTime.system(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )
    );

    const word fieldName = setSurfaceDistanceFieldDict.lookup("fieldName");

    Info<< "Reading field " << fieldName << "\n" << endl;
    volScalarField f
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    triSurfaceDistanceFunction surfDistFcn
    (
        mesh,
        setSurfaceDistanceFieldDict
    );

    Info<< "Evaluating field " << fieldName << "\n" << endl;
    surfDistFcn.evaluate(f);

    f.write();
}

