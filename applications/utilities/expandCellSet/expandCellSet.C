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
    expandCellSet

Description
    Expand a cell set by adding neighbours cells of cells in the set.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "cellSet.H"
#include "DynamicList.H"
#include "ListOps.H"
#include "autoPtr.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<cellSet> addNeighboursToSet
(
    const cellSet& baseSet,
    const fvMesh& mesh
)
{
    word expandedCellSetName = baseSet.name() + "Expanded";

    DynamicList<label> neighbourCells(6*baseSet.size());
    const labelListList& cellCells = mesh.cellCells();

    forAllConstIter(labelHashSet, baseSet, iter)
    {
        label cellI = iter.key();
        forAll(cellCells[cellI], j)
        {
            label neighbourI = cellCells[cellI][j];
            neighbourCells.append(neighbourI);
        }
    }

    labelList uniqueIndices;
    uniqueOrder(neighbourCells.shrink(), uniqueIndices);
    labelList uniqueNeighbours(neighbourCells, uniqueIndices);

    autoPtr<cellSet> pExpandedCellSet
    (
        new cellSet(mesh, expandedCellSetName, uniqueNeighbours)
    );

    pExpandedCellSet->addSet(baseSet);

    return pExpandedCellSet;
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "enlarge a cell set by adding all cell neighbours"
    );
    #include "addRegionOption.H"
    argList::addOption
    (
        "set",
        "cellSet",
        "cell set that should be expanded"
    );
    argList::addOption
    (
        "layers",
        "N",
        "number of neighbour layers that are added"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const bool setOpt   = args.optionFound("set");
    const bool layerOpt = args.optionFound("layers");

    word cellSetName;
    if (setOpt)
    {
        cellSetName = args.optionRead<word>("set");
    }
    else
    {
        cellSetName = "lumen";
    }

    label nLayers = 1;

    if (layerOpt)
    {
        nLayers = args.optionRead<label>("layers");
    }

    cellSet selectedSet(mesh, cellSetName);

    cellSet& previousSet = selectedSet;

    word expandedSetName;

    for(int i=0; i < nLayers; i++)
    {
        autoPtr<cellSet> pExpandedSet = addNeighboursToSet(previousSet, mesh);

        if (i == nLayers - 1)
        {
            pExpandedSet->write();
        }
        previousSet = pExpandedSet();
        expandedSetName = pExpandedSet->name();
    }

    Info<< "Wrote expanded cell set " << expandedSetName << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
