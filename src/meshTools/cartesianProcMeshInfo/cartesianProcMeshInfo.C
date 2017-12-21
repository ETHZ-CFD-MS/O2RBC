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

\*---------------------------------------------------------------------------*/

#include "cartesianProcMeshInfo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cartesianProcMeshInfo::cartesianProcMeshInfo(const polyMesh& mesh)
:
    mesh_(mesh),
    procBBoxes_(Pstream::nProcs())
{
    // Communicate bounding boxes
    if (Pstream::parRun())
    {
        // Get bounding box without reducing
        boundBox bBox(mesh.points(), false);

        List<point> minBoundList(Pstream::nProcs(), point::zero);
        List<point> maxBoundList(Pstream::nProcs(), point::zero);

        minBoundList[Pstream::myProcNo()] = bBox.min();
        maxBoundList[Pstream::myProcNo()] = bBox.max();

        reduce(minBoundList, sumOp<List<point> >());
        reduce(maxBoundList, sumOp<List<point> >());

        forAll(minBoundList, i)
        {
            procBBoxes_[i] = boundBox(minBoundList[i], maxBoundList[i]);
        }
    }
    else
    {
        procBBoxes_[0] = mesh.bounds();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

label
Foam::cartesianProcMeshInfo::findProcNo(const point& p) const
{
    label procI = -1;
    forAll(procBBoxes_, i)
    {
        if (procBBoxes_[i].contains(p))
        {
            return i;
        }
    }

    return procI;
}

Foam::labelList 
Foam::cartesianProcMeshInfo::findOverlappingProcIDs
(
    const boundBox& bBox
) const
{
    labelList procIDs;
    forAll(procBBoxes_, i)
    {
        if (procBBoxes_[i].overlaps(bBox))
        {
            procIDs.append(i);
        }
    }

    return procIDs;
}

// ************************************************************************* //
