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

#include "disconnectedZoneMesh.H"

#include "GeometricField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(disconnectedZoneMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::disconnectedZoneMesh::disconnectedZoneMesh(const IOobject& io)
:
    regionAddRemoveFvMesh(io)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::disconnectedZoneMesh::~disconnectedZoneMesh()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label
Foam::disconnectedZoneMesh::nZoneCells(const label zoneI) const
{
    return cellZones()[zoneI].size();
}


Foam::labelPairList
Foam::disconnectedZoneMesh::zoneBoundaryFaces(const label zoneI) const
{
    labelPairList boundaryFaces;
    const Map<label>& zoneMap = cellZones().zoneMap();

    forAll(boundary(), patchI)
    {
        const fvPatch& cPatch = boundary()[patchI];
        forAll(cPatch, faceI)
        {
            label ownerI = faceOwner()[cPatch.start() + faceI];
            if (zoneMap[ownerI] == zoneI)
            {
                boundaryFaces.append(labelPair(patchI, faceI));
            }
        }
    }

    return boundaryFaces;
}


Foam::boundBox
Foam::disconnectedZoneMesh::zoneBounds(const label zoneI) const
{
    return boundBox(pointField(points(), pointZones()[zoneI]));
}


bool
Foam::disconnectedZoneMesh::hasNonEmptyZone(const label zoneI) const
{
    if
    (
        0 <= zoneI && zoneI < pointZones().size()
     && 0 <= zoneI && zoneI < faceZones() .size()
     && 0 <= zoneI && zoneI < cellZones() .size()
    )
    {
        if
        (
            pointZones()[zoneI].size() > 0
         && faceZones ()[zoneI].size() > 0
         && cellZones ()[zoneI].size() > 0
        )
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    return false;
}


void
Foam::disconnectedZoneMesh::movePoints
(
    const pointField& pf,
    const label zoneI
)
{
    const label nZonePoints = cellZones()[zoneI].size();
    if (pf.size() != nZonePoints)
    {
        FatalErrorIn
        (
            "Foam::disconnectedZoneMesh::movePoints(const pointField&, "
            "const label)"
        ) << "The number of points in pf (" << pf.size() 
          << ") is not equal" << nl
          << "to the number of points in the cell zone " << zoneI
          << " (" << nZonePoints << ")."
          << abort(FatalError);
    }

    pointField transformedPts(points());

    UIndirectList<point>(transformedPts, cellZones()[zoneI]) = pf;

    fvMesh::movePoints(transformedPts);
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
