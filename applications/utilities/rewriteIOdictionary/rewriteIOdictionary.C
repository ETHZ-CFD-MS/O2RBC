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
    rewriteIOdictionary

Description
    Rewrite an IOdictionary to a file. This is useful for files with a lot
    of white spaces, such as those produced by PyFoam.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Rewrites an IOdictionary in the constant directory"
    );
    #include "addRegionOption.H"
    argList::addOption
    (
        "fileName",
        "file",
        "name of the file to rewrite"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    word fileName(args.optionRead<word>("fileName"));

    IOdictionary dict
    (
        IOobject
        (
            fileName,
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false // do not register in database
        )
    );
    Info<< "Finished reading dictionary" << nl;

    dict.regIOobject::write();

    Info<< "Wrote dictionary" << nl
        << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

