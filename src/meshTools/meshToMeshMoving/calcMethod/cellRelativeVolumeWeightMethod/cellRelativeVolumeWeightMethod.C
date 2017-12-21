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

\*---------------------------------------------------------------------------*/

#include "cellRelativeVolumeWeightMethod.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "addToRunTimeSelectionTable.H"

#include "tetOverlapVolume.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellRelativeVolumeWeightMethod, 1);
    addToRunTimeSelectionTable
    (
        myMeshToMeshMethod,
        cellRelativeVolumeWeightMethod,
        components
    );
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


void Foam::cellRelativeVolumeWeightMethod::calculateAddressing
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
)
{
    label srcCellI = srcSeedI;
    label tgtCellI = tgtSeedI;

    List<DynamicList<label> > srcToTgtAddr(src_.nCells());
    List<DynamicList<scalar> > srcToTgtWght(src_.nCells());

    List<DynamicList<label> > tgtToSrcAddr(tgt_.nCells());
    List<DynamicList<scalar> > tgtToSrcWght(tgt_.nCells());

    // list of tgt cell neighbour cells
    DynamicList<label> nbrTgtCells(10);

    // list of tgt cells currently visited for srcCellI to avoid multiple hits
    DynamicList<label> visitedTgtCells(10);

    // list to keep track of tgt cells used to seed src cells
    labelList seedCells(src_.nCells(), -1);
    seedCells[srcCellI] = tgtCellI;

    const scalarField& srcVol = src_.cellVolumes();
    const scalarField& tgtVol = tgt_.cellVolumes();

    label nCalls = 0;
    label nCallsZero = 0;

    // loop over candidate source cells
    do
    {
        nbrTgtCells.clear();
        visitedTgtCells.clear();

        // append initial target cell and neighbours
        nbrTgtCells.append(tgtCellI);
        appendNbrCells(tgtCellI, tgt_, visitedTgtCells, nbrTgtCells);

        // loop over candidate target cells
        do
        {
            tgtCellI = nbrTgtCells.remove();
            visitedTgtCells.append(tgtCellI);

            // check whether the (Cartesian) target cell intersects the current
            // source cell
            if (intersect(srcCellI, tgtCellI))
            {
                scalar vol = interVol(srcCellI, tgtCellI);
                nCalls++;
                if (vol == 0)
                {
                    nCallsZero++;
                }

                // accumulate addressing and weights for valid intersection
                if (vol/srcVol[srcCellI] > tolerance_)
                {
                    // store src/tgt cell pair
                    srcToTgtAddr[srcCellI].append(tgtCellI);
                    
                    // divide intersection volume by the volume of the donor cell
                    srcToTgtWght[srcCellI].append(vol/srcVol[srcCellI]);

                    tgtToSrcAddr[tgtCellI].append(srcCellI);
                    
                    // divide intersection volume by the volume of the donor cell
                    tgtToSrcWght[tgtCellI].append(vol/tgtVol[tgtCellI]);

                    appendNbrCells(tgtCellI, tgt_, visitedTgtCells, nbrTgtCells);

                    // accumulate intersection volume
                    V_ += vol;
                }
            }
        }
        while (!nbrTgtCells.empty());

        mapFlag[srcCellI] = false;

        // find new source seed cell
        setNextCells
        (
            startSeedI,
            srcCellI,
            tgtCellI,
            srcCellIDs,
            mapFlag,
            visitedTgtCells,
            seedCells
        );
    }
    while (srcCellI != -1);

    // transfer addressing into persistent storage
    forAll(srcToTgtCellAddr, i)
    {
        srcToTgtCellAddr[i].transfer(srcToTgtAddr[i]);
        srcToTgtCellWght[i].transfer(srcToTgtWght[i]);
    }

    forAll(tgtToSrcCellAddr, i)
    {
        tgtToSrcCellAddr[i].transfer(tgtToSrcAddr[i]);
        tgtToSrcCellWght[i].transfer(tgtToSrcWght[i]);
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellRelativeVolumeWeightMethod::cellRelativeVolumeWeightMethod
(
    const polyMesh& src,
    const polyMesh& tgt
)
:
    myCellVolumeWeightMethod(src, tgt)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::cellRelativeVolumeWeightMethod::~cellRelativeVolumeWeightMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellRelativeVolumeWeightMethod::calculate
(
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght,
    labelListList& tgtToSrcAddr,
    scalarListList& tgtToSrcWght
)
{
    bool ok = initialise
    (
        srcToTgtAddr,
        srcToTgtWght,
        tgtToSrcAddr,
        tgtToSrcWght
    );

    if (!ok)
    {
        return;
    }

    // (potentially) participating source mesh cells
    const labelList srcCellIDs(maskCells());

    // list to keep track of whether src cell can be mapped
    boolList mapFlag(src_.nCells(), false);
    UIndirectList<bool>(mapFlag, srcCellIDs) = true;

    // find initial point in tgt mesh
    label srcSeedI = -1;
    label tgtSeedI = -1;
    label startSeedI = 0;

    bool startWalk =
        findInitialSeeds
        (
            srcCellIDs,
            mapFlag,
            startSeedI,
            srcSeedI,
            tgtSeedI
        );

    if (startWalk)
    {
        calculateAddressing
        (
            srcToTgtAddr,
            srcToTgtWght,
            tgtToSrcAddr,
            tgtToSrcWght,
            srcSeedI,
            tgtSeedI,
            srcCellIDs,
            mapFlag,
            startSeedI
        );
    }
    else
    {
        WarningInFunction
            << "Could not find an initial seed. Setting empty addressing "
            << "and weight lists." << endl;
        tgtToSrcAddr = labelListList(tgt_.nCells(), labelList());
        tgtToSrcWght = scalarListList(tgt_.nCells(), scalarList());
        srcToTgtAddr = labelListList(src_.nCells(), labelList());
        srcToTgtWght = scalarListList(src_.nCells(), scalarList());

        // if meshes are collocated, after inflating the source mesh bounding
        // box tgt mesh cells may be transferred, but may still not overlap
        // with the source mesh
        return;
    }
}


// ************************************************************************* //
