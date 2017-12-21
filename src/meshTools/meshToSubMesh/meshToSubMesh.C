/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "meshToSubMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

meshToSubMesh::meshToSubMesh(const fvMesh& mesh, const fvMesh& subMesh)
:
    mesh_(mesh),
    subMesh_(subMesh),
    cellIdSubMeshToMesh_(),
    faceIdSubMeshToMesh_(),
    boundaryIdSubMeshToMesh_()
{
    createLists();
}

meshToSubMesh::meshToSubMesh
(
    const meshToSubMesh& obj
)
:
    mesh_(obj.mesh_),
    subMesh_(obj.subMesh_),
    cellIdSubMeshToMesh_(obj.cellIdSubMeshToMesh_),
    faceIdSubMeshToMesh_(obj.faceIdSubMeshToMesh_),
    boundaryIdSubMeshToMesh_(obj.boundaryIdSubMeshToMesh_)
{}


void meshToSubMesh::createLists()
{
    // read list of indices:
    // i-th element in the list => "global" index of i-th cell in subMesh

    labelIOList cellAddressing = 
        IOobject
        (
            "cellRegionAddressing",
            subMesh_.facesInstance(),
            subMesh_.meshSubDir,
            subMesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
    
    cellIdSubMeshToMesh_ = cellAddressing;

    labelIOList faceAddressing = 
        IOobject
        (
            "faceRegionAddressing",
            subMesh_.facesInstance(),
            subMesh_.meshSubDir,
            subMesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

    // The index list in faceRegionAddressing takes orientation into account as
    // follows (see http://openfoamwiki.net/index.php/SplitMeshRegions)
    // faceRegionAddressing : for every face in this region the face in the original mesh 
    // + "turning index". For a face in the same orientation this is the original facelabel+1, 
    // for a turned face this is -facelabel-1
    forAll(faceAddressing, j)
    {
        faceAddressing[j] = abs(faceAddressing[j]) - 1;
    }
    faceIdSubMeshToMesh_ = faceAddressing;

    labelIOList boundaryAddressing = 
        IOobject
        (
            "boundaryRegionAddressing",
            subMesh_.facesInstance(),
            subMesh_.meshSubDir,
            subMesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
    boundaryIdSubMeshToMesh_ = boundaryAddressing;

}


void meshToSubMesh::write(Ostream& os) const
{
    os << " Copy object between mesh " << mesh_.name() 
       << " and mesh " << subMesh_.name() << endl;
}


// ************************************************************************* //
