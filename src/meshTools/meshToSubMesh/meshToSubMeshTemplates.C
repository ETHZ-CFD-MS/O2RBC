/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void meshToSubMesh::copyBoundaryFieldMeshToSubMesh
(
    const GeometricField<Type,fvsPatchField,surfaceMesh>& fFaces,
    GeometricField<Type,PatchField,GeoMesh>& fSub
) const
{
    // copy values on patches
    const fvPatchList& patches = subMesh_.boundary();

    forAll(patches, patchI)
    {
        const fvPatch& cPatch = patches[patchI];

        // ID of the boundary in the global mesh corresponding to cPatch
        const label meshBoundaryID = boundaryIdSubMeshToMesh_[patchI];
        
        // If current patch is a patch field of the global mesh,
        // access the boundary field of the global mesh. 
        if (meshBoundaryID >= 0)
        {
            const fvPatch& domainPatch = mesh_.boundary()[meshBoundaryID];

            forAll(cPatch, faceI)
            {
                // index of the face in subMesh
                label subMeshFaceIndex = cPatch.start() + faceI;
                // index of the face in the global mesh
                label meshFaceIndex = faceIdSubMeshToMesh_[subMeshFaceIndex];
                // local index of the face on the patch in the global mesh
                label meshBoundaryFaceIndex = meshFaceIndex - domainPatch.start();
                // access the corresponding boundary value of fFaces
                fSub.boundaryFieldRef()[patchI][faceI] = 
                        fFaces.boundaryField()[meshBoundaryID][meshBoundaryFaceIndex];
            }
        }
        // If current patch is a not a patch field of the main mesh,
        // it is part of an internal field of the main mesh.
        // Thus the internal field of fFaces is accessed.
        else
        {
            forAll(cPatch, faceI)
            {
                // index of the face in subMesh
                label subMeshFaceIndex = cPatch.start() + faceI;
                // index of the face in the global mesh
                label meshFaceIndex = faceIdSubMeshToMesh_[subMeshFaceIndex];
                // get the corresponding value from the internal field of fFaces
                fSub.boundaryFieldRef()[patchI][faceI] = fFaces[meshFaceIndex];
            }
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void meshToSubMesh::copyBoundaryFieldSubMeshToMesh
(
    GeometricField<Type,PatchField,GeoMesh>& f,
    const GeometricField<Type,PatchField,GeoMesh>& fSub
) const
{
    // copy values on patches
    const fvPatchList& patches = subMesh_.boundary();

    forAll(patches, patchI)
    {
        const fvPatch& cPatch = patches[patchI];

        // ID of the boundary in the global mesh corresponding to cPatch
        const label meshBoundaryID = boundaryIdSubMeshToMesh_[patchI];
        
        // If current patch is a patch field of the global mesh,
        // access the boundary field of the global mesh. 
        if (meshBoundaryID >= 0)
        {
            const fvPatch& domainPatch = mesh_.boundary()[meshBoundaryID];

            forAll(cPatch, faceI)
            {
                // index of the face in subMesh
                label subMeshFaceIndex = cPatch.start() + faceI;
                // index of the face in the global mesh
                label meshFaceIndex = faceIdSubMeshToMesh_[subMeshFaceIndex];
                // local index of the face on the patch in the global mesh
                label meshBoundaryFaceIndex = meshFaceIndex - domainPatch.start();
                // access the corresponding boundary value of fSub
                f.boundaryFieldRef()[meshBoundaryID][meshBoundaryFaceIndex] =
                        fSub.boundaryField()[patchI][faceI];
            }
        }
    }
}


template<class Type>
void meshToSubMesh::copySubMeshToMeshCells
(
    Field<Type>& f,
    const Field<Type>& fSub
) const
{
    // copy cell-centered values
    forAll(fSub, i)
    {
        f[cellIdSubMeshToMesh_[i]] = fSub[i];
    }
}


template<class Type>
void meshToSubMesh::copyMeshToSubMeshCells
(
    const Field<Type>& f,
    Field<Type>& fSub
) const
{
    // copy cell-centered values
    forAll(fSub, i)
    {
        fSub[i] = f[cellIdSubMeshToMesh_[i]];
    }
}


template<class Type>
void meshToSubMesh::copySubMeshToMeshFaces
(
    Field<Type>& f,
    const Field<Type>& fSub
) const
{
    // copy face-centered values
    forAll(fSub, i)
    {
        f[faceIdSubMeshToMesh_[i]] = fSub[i];
    }
}


template<class Type>
void meshToSubMesh::copyMeshToSubMeshFaces
(
    const Field<Type>& f,
    Field<Type>& fSub
) const
{
    // copy face-centered values
    forAll(fSub, i)
    {
        fSub[i] = f[faceIdSubMeshToMesh_[i]];
    }
}


template<class Type>
void meshToSubMesh::copySubMeshToMesh
(
    GeometricField<Type,fvPatchField,volMesh>& f,
    const GeometricField<Type,fvPatchField,volMesh>& fSub
) const
{
    copySubMeshToMeshCells(f.internalFieldRef(), fSub.internalField());

    // copy values on patches
    copyBoundaryFieldSubMeshToMesh(f, fSub);
}


template<class Type>
void meshToSubMesh::copyMeshToSubMesh
(
    const GeometricField<Type,fvPatchField,volMesh>& f,
    GeometricField<Type,fvPatchField,volMesh>& fSub
) const
{
    copyMeshToSubMeshCells(f.internalField(), fSub.internalFieldRef());

    // copy values on patches
    surfaceScalarField fFaces = linearInterpolate(f);
    copyBoundaryFieldMeshToSubMesh(fFaces, fSub);
}


template<class Type>
void meshToSubMesh::copySubMeshToMesh
(
    GeometricField<Type,fvsPatchField,surfaceMesh>& f,
    const GeometricField<Type,fvsPatchField,surfaceMesh>& fSub
) const
{
    copySubMeshToMeshFaces(f.primitiveFieldRef(), fSub.internalField());

    // copy values on patches
    copyBoundaryFieldSubMeshToMesh(f, fSub);
}


template<class Type>
void meshToSubMesh::copyMeshToSubMesh
(
    const GeometricField<Type,fvsPatchField,surfaceMesh>& f,
    GeometricField<Type,fvsPatchField,surfaceMesh>& fSub
) const
{
    copyMeshToSubMeshFaces(f.internalField(), fSub.primitiveFieldRef());

    // copy values on patches
    copyBoundaryFieldMeshToSubMesh(f, fSub);
}

