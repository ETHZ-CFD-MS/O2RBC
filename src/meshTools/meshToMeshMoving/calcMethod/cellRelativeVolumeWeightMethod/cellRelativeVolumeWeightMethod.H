/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
    Foam::cellRelativeVolumeWeightMethod

Description
    Cell-volume-weighted mesh-to-mesh interpolation class

    Volume conservative. The interpolation weights are defined by the volume
    ratio of the source cell to the target cell (or vice versa).

SourceFiles
    cellRelativeVolumeWeightMethod.C

\*---------------------------------------------------------------------------*/

#ifndef cellRelativeVolumeWeightMethod_H
#define cellRelativeVolumeWeightMethod_H

#include "myCellVolumeWeightMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cellRelativeVolumeWeightMethod Declaration
\*---------------------------------------------------------------------------*/

class cellRelativeVolumeWeightMethod
:
    public myCellVolumeWeightMethod
{
protected:

    // Protected Member Functions

        //- Calculate the mesh-to-mesh addressing and weights
        void calculateAddressing
        (
            labelListList& srcToTgtCellAddr,
            scalarListList& srcToTgtCellWght,
            labelListList& tgtToSrcCellAddr,
            scalarListList& tgtToSrcCellWght,
            const label srcSeedI,
            const label tgtSeedI,
            const labelList& srcCellIDs,
            boolList& mapFlag,
            label& startSeedI
        );

        //- Disallow default bitwise copy construct
        cellRelativeVolumeWeightMethod(const cellRelativeVolumeWeightMethod&);

        //- Disallow default bitwise assignment
        void operator=(const cellRelativeVolumeWeightMethod&);


public:

    //- Run-time type information
    TypeName("cellRelativeVolumeWeight");

    //- Construct from source and target meshes
    cellRelativeVolumeWeightMethod(const polyMesh& src, const polyMesh& tgt);

    //- Destructor
    virtual ~cellRelativeVolumeWeightMethod();


    // Member Functions

        // Evaluate

            //- Calculate addressing and weights
            virtual void calculate
            (
                labelListList& srcToTgtAddr,
                scalarListList& srcToTgtWght,
                labelListList& tgtToTgtAddr,
                scalarListList& tgtToTgtWght
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
