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

\*---------------------------------------------------------------------------*/

#include "regionAddRemoveFvMesh.H"
#include "Time.H"
#include "polyTopoChanger.H"
#include "mapPolyMesh.H"
#include "polyAddPoint.H"
#include "polyAddCell.H"
#include "polyAddFace.H"
#include "polyRemovePoint.H"
#include "polyRemoveCell.H"
#include "polyRemoveFace.H"
#include "polyModifyFace.H"

#include "GeometricField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionAddRemoveFvMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::polyTopoChange&
Foam::regionAddRemoveFvMesh::meshModifier() const
{
    if (meshModPtr_.empty())
    {
        meshModPtr_.reset
        (
            new polyTopoChange(*this)
        );
    }
    return meshModPtr_();
}


Foam::label Foam::regionAddRemoveFvMesh::patchIndex(const polyPatch& p)
{
    // Find the patch name on the list.  If the patch is already there
    // and patch types match, return index
    const word& pType = p.type();
    const word& pName = p.name();

    bool nameFound = false;

    forAll(patchNames_, patchI)
    {
        if (patchNames_[patchI] == pName)
        {
            if (word(patchDicts_[patchI]["type"]) == pType)
            {
                // Found name and types match
                return patchI;
            }
            else
            {
                // Found the name, but type is different
                nameFound = true;
            }
        }
    }

    // Patch not found.  Append to the list
    {
        OStringStream os;
        p.write(os);
        patchDicts_.append(dictionary(IStringStream(os.str())()));
    }

    if (nameFound)
    {
        // Duplicate name is not allowed.  Create a composite name from the
        // patch name and case name
        const word& caseName = p.boundaryMesh().mesh().time().caseName();

        patchNames_.append(pName + "_" + caseName);

        Info<< "label patchIndex(const polyPatch& p) : "
            << "Patch " << p.index() << " named "
            << pName << " in mesh " << caseName
            << " already exists, but patch types "
            << " do not match.\nCreating a composite name as "
            << patchNames_.last() << endl;
    }
    else
    {
        patchNames_.append(pName);
    }

    return patchNames_.size() - 1;
}


Foam::label Foam::regionAddRemoveFvMesh::zoneIndex
(
    DynamicList<word>& names,
    const word& curName
)
{
    forAll(names, zoneI)
    {
        if (names[zoneI] == curName)
        {
            return zoneI;
        }
    }

    // Not found.  Add new name to the list
    names.append(curName);

    return names.size() - 1;
}


void Foam::regionAddRemoveFvMesh::changePatchID
(
    const label faceID,
    const label patchID
)
{
    const label zoneID = this->faceZones().whichZone(faceID);

    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = this->faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
    }

    polyTopoChange& meshMod = meshModifier();

    meshMod.setAction
    (
        polyModifyFace
        (
            this->faces()[faceID],              // face
            faceID,                             // face ID
            this->faceOwner()[faceID],          // owner
            -1,                                 // neighbour
            false,                              // flip flux
            patchID,                            // patch ID
            false,                              // remove from zone
            zoneID,                             // zone ID
            zoneFlip                            // zone flip
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionAddRemoveFvMesh::regionAddRemoveFvMesh(const IOobject& io)
:
    topoChangerFvMesh(io),
    meshModPtr_(NULL),
    patchNames_(2*boundaryMesh().size()),
    patchDicts_(2*boundaryMesh().size()),
    pointZoneNames_(),
    faceZoneNames_(),
    cellZoneNames_(),
    meshMap_()
{
    // Insert the original patches into the list
    wordList curPatchNames = boundaryMesh().names();

    forAll(boundaryMesh(), patchI)
    {
        patchNames_.append(boundaryMesh()[patchI].name());

        OStringStream os;
        boundaryMesh()[patchI].write(os);
        patchDicts_.append(dictionary(IStringStream(os.str())()));
    }

    // Insert point, face and cell zones into the list

    // Point zones
    wordList curPointZoneNames = pointZones().names();
    if (curPointZoneNames.size())
    {
        pointZoneNames_.setCapacity(2*curPointZoneNames.size());
    }

    forAll(curPointZoneNames, zoneI)
    {
        pointZoneNames_.append(curPointZoneNames[zoneI]);
    }

    // Face zones
    wordList curFaceZoneNames = faceZones().names();

    if (curFaceZoneNames.size())
    {
        faceZoneNames_.setCapacity(2*curFaceZoneNames.size());
    }
    forAll(curFaceZoneNames, zoneI)
    {
        faceZoneNames_.append(curFaceZoneNames[zoneI]);
    }

    // Cell zones
    wordList curCellZoneNames = cellZones().names();

    if (curCellZoneNames.size())
    {
        cellZoneNames_.setCapacity(2*curCellZoneNames.size());
    }
    forAll(curCellZoneNames, zoneI)
    {
        cellZoneNames_.append(curCellZoneNames[zoneI]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionAddRemoveFvMesh::~regionAddRemoveFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionAddRemoveFvMesh::addMesh
(
    const polyMesh& m
)
{
    polyTopoChange& meshMod = meshModifier();

    // Add all the points, faces and cells of the new mesh

    // Add points

    label zoneID = -1;

    const pointField& p = m.points();
    labelList renumberPoints(p.size());

    const pointZoneMesh& pz = m.pointZones();
    labelList pointZoneIndices(pz.size());

    forAll(pz, zoneI)
    {
        pointZoneIndices[zoneI] = zoneIndex(pointZoneNames_, pz[zoneI].name());
    }

    forAll(p, pointI)
    {
        // Grab zone ID.  If a point is not in a zone, it will return -1
        zoneID = pz.whichZone(pointI);

        if (zoneID >= 0)
        {
            // Translate zone ID into the new index
            zoneID = pointZoneIndices[zoneID];
        }

        renumberPoints[pointI] =
            meshMod.setAction
            (
                polyAddPoint
                (
                    p[pointI],            // Point to add
                    -1,                   // Master point (straight addition)
                    zoneID,               // Zone for point
                    pointI < m.nPoints()  // Is in cell?
                )
            );
    }

    // Add cells

    const cellList& c = m.cells();
    labelList renumberCells(c.size());

    const cellZoneMesh& cz = m.cellZones();
    labelList cellZoneIndices(cz.size());

    forAll(cz, zoneI)
    {
        cellZoneIndices[zoneI] = zoneIndex(cellZoneNames_, cz[zoneI].name());
    }

    forAll(c, cellI)
    {
        // Grab zone ID.  If a cell is not in a zone, it will return -1
        zoneID = cz.whichZone(cellI);

        if (zoneID >= 0)
        {
            // Translate zone ID into the new index
            zoneID = cellZoneIndices[zoneID];
        }

        renumberCells[cellI] =
            meshMod.setAction
            (
                polyAddCell
                (
                    -1,                   // Master point
                    -1,                   // Master edge
                    -1,                   // Master face
                    -1,                   // Master cell
                    zoneID                // Zone for cell
                )
            );
    }

    // Add faces
    const polyBoundaryMesh& bm = m.boundaryMesh();

    // Gather the patch indices
    labelList patchIndices(bm.size());

    forAll(patchIndices, patchI)
    {
        patchIndices[patchI] = patchIndex(bm[patchI]);
    }

    // Temporary: update number of allowable patches. This should be
    // determined at the top - before adding anything.
    meshMod.setNumPatches(patchNames_.size());

    const faceZoneMesh& fz = m.faceZones();
    labelList faceZoneIndices(fz.size());

    forAll(fz, zoneI)
    {
        faceZoneIndices[zoneI] = zoneIndex(faceZoneNames_, fz[zoneI].name());
    }

    const faceList& f = m.faces();
    labelList renumberFaces(f.size());

    const labelList& own = m.faceOwner();
    const labelList& nei = m.faceNeighbour();

    label newOwn, newNei, newPatch, newZone;
    bool newZoneFlip;

    forAll(f, faceI)
    {
        const face& curFace = f[faceI];

        face newFace(curFace.size());

        forAll(curFace, pointI)
        {
            newFace[pointI] = renumberPoints[curFace[pointI]];
        }

        if (debug)
        {
            // Check that the face is valid
            if (min(newFace) < 0)
            {
                FatalErrorIn("void regionAddRemoveFvMesh::addMesh(const polyMesh&)")
                    << "Error in point mapping for face " << faceI
                    << ".  Old face: " << curFace << " New face: " << newFace
                    << abort(FatalError);
            }
        }

        if (faceI < m.nInternalFaces() || faceI >= m.nFaces())
        {
            newPatch = -1;
        }
        else
        {
            newPatch = patchIndices[bm.whichPatch(faceI)];
        }

        newOwn = own[faceI];
        if (newOwn > -1) newOwn = renumberCells[newOwn];

        if (newPatch > -1)
        {
            newNei = -1;
        }
        else
        {
            newNei = nei[faceI];
            newNei = renumberCells[newNei];
        }


        newZone = fz.whichZone(faceI);
        newZoneFlip = false;

        if (newZone >= 0)
        {
            newZoneFlip = fz[newZone].flipMap()[fz[newZone].whichFace(faceI)];

            // Grab the new zone
            newZone = faceZoneIndices[newZone];
        }

        renumberFaces[faceI] =
            meshMod.setAction
            (
                polyAddFace
                (
                    newFace,
                    newOwn,
                    newNei,
                    -1,
                    -1,
                    -1,
                    false,
                    newPatch,
                    newZone,
                    newZoneFlip
                )
            );
    }

}


void Foam::regionAddRemoveFvMesh::removeMesh
(
    const word& cellZoneName,
    const word& pointZoneName,
    const word& faceZoneName
)
{
    polyTopoChange& meshMod = meshModifier();

    // Check that there is a pointZone, faceZone and cellZone
    // with one of the zone names
    if (cellZones().findZoneID(cellZoneName) == -1 
     || pointZones().findZoneID(pointZoneName) == -1
     || faceZones().findZoneID(faceZoneName) == -1)
    {
        FatalErrorIn
        (
            "void regionAddRemoveFvMesh::removeMesh(const word&)"
        ) << "Cell zone, point zone or face zone with name "
          << cellZoneName << ", " << pointZoneName << ", "
          << faceZoneName << " not found." << nl
          << abort(FatalError);
    }

    // Remove points
    label zoneID = pointZones().findZoneID(pointZoneName);

    const pointZone& pz = pointZones()[zoneID];
    forAll(pz, i)
    {
        label pI = pz[i];

        meshMod.setAction(polyRemovePoint(pI));
    }

    // Remove cells
    zoneID = cellZones().findZoneID(cellZoneName);

    const cellZone& cz = cellZones()[zoneID];
    forAll(cz, i)
    {
        label cI = cz[i];

        meshMod.setAction(polyRemoveCell(cI));
    }

    // Remove faces
    zoneID = faceZones().findZoneID(faceZoneName);

    const faceZone& fz = faceZones()[zoneID];
    forAll(fz, i)
    {
        label fI = fz[i];

        meshMod.setAction(polyRemoveFace(fI));
    }
}


void Foam::regionAddRemoveFvMesh::merge()
{
    if (debug)
    {
        Info<< "regionAddRemoveFvMesh::merge() :" << nl
            << "    patch names: " << patchNames_ << nl
            << "    patch dicts: " << patchDicts_ << nl
            << "    point zone names: " << pointZoneNames_ << nl
            << "    face zone names: " << faceZoneNames_ << nl
            << "    cell zone names: " << cellZoneNames_ << endl;
    }

    // Add the patches if necessary
    if (patchNames_.size() != boundaryMesh().size())
    {
        Info<< "Copying old patches" << endl;

        List<polyPatch*> newPatches(patchNames_.size());

        const polyBoundaryMesh& oldPatches = boundaryMesh();

        // Note.  Re-using counter in two for loops
        label patchI = 0;

        for (patchI = 0; patchI < oldPatches.size(); patchI++)
        {
            newPatches[patchI] = oldPatches[patchI].clone(oldPatches).ptr();
        }

        Info<< "Adding new patches. " << endl;
        WarningIn
        (
            "Foam::regionAddRemoveFvMesh::merge()"
        ) << "Patches need to be added. Check whether this is "
          << "the desired behaviour." << nl
          << "patchNames_ = " << patchNames_
          << "original patch names: " << boundaryMesh().names() << nl << endl;

        label endOfLastPatch =
            oldPatches[patchI - 1].start() + oldPatches[patchI - 1].size();

        for (; patchI < patchNames_.size(); patchI++)
        {
            // Add a patch
            dictionary dict(patchDicts_[patchI]);
            dict.set("nFaces", 0);
            dict.set("startFace", endOfLastPatch);

            newPatches[patchI] =
            (
                polyPatch::New
                (
                    patchNames_[patchI],
                    dict,
                    patchI,
                    oldPatches
                ).ptr()
            );
        }

        removeBoundary();
        addPatches(newPatches);
    }

    // Add the zones if necessary
    if (pointZoneNames_.size() > pointZones().size())
    {
        label nZones = pointZones().size();

        pointZones().setSize(pointZoneNames_.size());

        for (label zoneI = nZones; zoneI < pointZoneNames_.size(); zoneI++)
        {
            pointZones().set
            (
                zoneI,
                new pointZone
                (
                    pointZoneNames_[zoneI],
                    labelList(),
                    zoneI,
                    pointZones()
                )
            );
        }
    }
    if (cellZoneNames_.size() > cellZones().size())
    {
        label nZones = cellZones().size();

        cellZones().setSize(cellZoneNames_.size());

        for (label zoneI = nZones; zoneI < cellZoneNames_.size(); zoneI++)
        {
            cellZones().set
            (
                zoneI,
                new cellZone
                (
                    cellZoneNames_[zoneI],
                    labelList(),
                    zoneI,
                    cellZones()
                )
            );
        }
    }
    if (faceZoneNames_.size() > faceZones().size())
    {
        label nZones = faceZones().size();

        faceZones().setSize(faceZoneNames_.size());

        for (label zoneI = nZones; zoneI < faceZoneNames_.size(); zoneI++)
        {
            faceZones().set
            (
                zoneI,
                new faceZone
                (
                    faceZoneNames_[zoneI],
                    labelList(),
                    boolList(),
                    zoneI,
                    faceZones()
                )
            );
        }
    }
}


void
Foam::regionAddRemoveFvMesh::filterPatches()
{
    polyTopoChange& meshMod = meshModifier();

    const polyBoundaryMesh& patches = this->boundaryMesh();
    label destPatchI = patches.findPatchID("RBCBoundary");

    // remove the empty patches
    const polyPatch& uniquePatch = patches[destPatchI];
    List<polyPatch*> onePatchList(1);
    onePatchList[0] = 
        uniquePatch.clone
        (
            patches,
            1,
            uniquePatch.size(),
            uniquePatch.start()
        ).ptr();

    this->removeBoundary();
    this->addPatches(onePatchList);

    // tell the polyTopoChange object that there is only one patch left
    meshMod.setNumPatches(1);

    Info<< "Finished filtering patches..." << endl;
}


bool
Foam::regionAddRemoveFvMesh::update()
{
    polyTopoChange& meshMod = meshModifier();

    // Change mesh. No inflation
    meshMap_ = meshMod.changeMesh
               (
                    *this, 
                    /*inflate*/ false, 
                    /*syncParallel*/ false
               );

    // Update mesh (this update the fields as well)
    if (debug)
    {
        Pout<< "Foam::regionAddRemoveFvMesh::update(): "
            << "Calling updateMesh(const mapPolyMesh&)" << endl;
    }

    this->updateMesh(meshMap_);

    if (debug)
    {
        Pout<< "Foam::regionAddRemoveFvMesh::update(): "
            << "Finished updateMesh(const mapPolyMesh&)" << endl;
    }

    // Clear topo changes for the next operation
    meshModPtr_.clear();

    return true;
}


// ************************************************************************* //

