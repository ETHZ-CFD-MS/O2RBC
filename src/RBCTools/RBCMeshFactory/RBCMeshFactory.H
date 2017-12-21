/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

Class
    Foam::RBCMeshFactory

Description
    "Factory" for the creation of a polyMesh for a red blood cell. The created
    polyMesh has a unique pointZone, faceZone and cellZone with names given
    by a prefix and an index.

SourceFiles
    RBCMeshFactory.C

\*---------------------------------------------------------------------------*/

#ifndef RBCMeshFactory_H
#define RBCMeshFactory_H

#include "fvCFD.H"

#include "deformableBodyGeometricState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class RBCMeshFactory Declaration
\*---------------------------------------------------------------------------*/

class RBCMeshFactory
{
    // Static private data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        RBCMeshFactory(const RBCMeshFactory&);

        //- Disallow default bitwise assignment
        void operator=(const RBCMeshFactory&);


public:

    // Static data

        //- Prefix for zone names
        static const word RBCPrefix;

    // Static member functions

        //- Create a polyMesh from a directory and an index number
        static autoPtr<polyMesh> createRBCMesh
        (
            const word&,
            const Time&,
            const label
        );

        //- Return a point zone name from an index number
        static word pointZoneName(label);

        //- Return a cell zone name from an index number
        static word cellZoneName(label);

        //- Return a face zone name from an index number
        static word faceZoneName(label);

        //- Return the geometric state of the mesh passed as argument
        static autoPtr<deformableBodyGeometricState> geometricState
        (
            const polyMesh& mesh
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
