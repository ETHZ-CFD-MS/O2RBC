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
    Foam::cartesianProcMeshInfo

Description
    Gathers information about the processors meshes. Works under the assumption
    that the volumes spanned by processor meshes are equal to their bounding 
    box, which is satisfied by Cartesian meshes.

SourceFiles
    cartesianProcMeshInfo.C

\*---------------------------------------------------------------------------*/

#ifndef cartesianProcMeshInfo_H
#define cartesianProcMeshInfo_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cartesianProcMeshInfo Declaration
\*---------------------------------------------------------------------------*/

class cartesianProcMeshInfo
{
    // Private data

        //- Reference to the mesh
        const polyMesh& mesh_;

        //- Bounding boxes of all processor meshes
        List<boundBox> procBBoxes_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cartesianProcMeshInfo(const cartesianProcMeshInfo&);

        //- Disallow default bitwise assignment
        void operator=(const cartesianProcMeshInfo&);


public:

    // Constructors

        //- Construct from components
        cartesianProcMeshInfo(const polyMesh&);

    // Member Functions

        // Access
            const polyMesh& mesh() const
            {
                return mesh_;
            }

            const List<boundBox>& procBBoxes() const
            {
                return procBBoxes_;
            }

            //- Returns the index of the processor that contains the argument.
            //  If the point is outside all processor meshes, return -1.
            label findProcNo(const point&) const;

            //- Return a list with the processor indices whose bounding box
            //  overlaps the bounding box given as argument.
            labelList findOverlappingProcIDs(const boundBox&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
