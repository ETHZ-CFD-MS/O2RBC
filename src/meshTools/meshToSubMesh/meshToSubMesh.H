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

Class
    Foam::meshToSubMesh

Description
    This class handles the copy of fields from a given mesh to a submesh.

SourceFiles
    meshToSubMesh.C

\*---------------------------------------------------------------------------*/

#ifndef meshToSubMesh_H
#define meshToSubMesh_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class meshToSubMesh
{
    // Private data

        //- Mesh object
        const fvMesh& mesh_;

        //- Mesh that consists of cells of mesh_.
        //  Typically produced by splitMeshRegions.
        const fvMesh& subMesh_;

        //- Addressing list for cells
        labelList cellIdSubMeshToMesh_;

        //- Addressing list for faces
        labelList faceIdSubMeshToMesh_;

        //- Addressing list for the boundary
        labelList boundaryIdSubMeshToMesh_;

    // Private member functions

        //- Create the addressing lists used for field copy
        void createLists();

        //- Copy patch field from mesh to subMesh
        template<class Type, template<class> class PatchField, class GeoMesh>
        void copyBoundaryFieldMeshToSubMesh
        (
            const GeometricField<Type,fvsPatchField,surfaceMesh>& fFaces,
            GeometricField<Type,PatchField,GeoMesh>& fSub
        ) const;

        //- Copy patch field from subMesh to mesh
        template<class Type, template<class> class PatchField, class GeoMesh>
        void copyBoundaryFieldSubMeshToMesh
        (
            GeometricField<Type,PatchField,GeoMesh>& f,
            const GeometricField<Type,PatchField,GeoMesh>& fSub
        ) const;

public:
    // Constructors

        //- Construct from attributes
        meshToSubMesh
        (
            const fvMesh& mesh, 
            const fvMesh& subMesh
        );

        //- Copy constructor
        meshToSubMesh
        (
            const meshToSubMesh&
        );

    // Member functions

        // Access
            
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            const fvMesh& subMesh() const
            {
                return subMesh_;
            }

        // Copy functions

            //- Copy a field defined on cells from the subMesh to the mesh
            template<class Type>
            void copySubMeshToMeshCells
            (
                Field<Type>& f,
                const Field<Type>& fSub
            ) const;

            //- Copy a field defined on cells defined on the mesh to the subMesh
            template<class Type>
            void copyMeshToSubMeshCells
            (
                const Field<Type>& f,
                Field<Type>& fSub
            ) const;

            //- Copy a field defined on faces from the subMesh to the mesh
            template<class Type>
            void copySubMeshToMeshFaces
            (
                Field<Type>& f,
                const Field<Type>& fSub
            ) const;

            //- Copy a field defined on faces defined on the mesh to the subMesh
            template<class Type>
            void copyMeshToSubMeshFaces
            (
                const Field<Type>& f,
                Field<Type>& fSub
            ) const;

            //- Copy a volumetric field from the subMesh to the mesh
            template<class Type>
            void copySubMeshToMesh
            (
                GeometricField<Type,fvPatchField,volMesh>& f,
                const GeometricField<Type,fvPatchField,volMesh>& fSub
            ) const;

            //- Copy a volumetric field defined on the mesh to the subMesh
            template<class Type>
            void copyMeshToSubMesh
            (
                const GeometricField<Type,fvPatchField,volMesh>& f,
                GeometricField<Type,fvPatchField,volMesh>& fSub
            ) const;

            //- Copy a surface field from the subMesh to the mesh
            template<class Type>
            void copySubMeshToMesh
            (
                GeometricField<Type,fvsPatchField,surfaceMesh>& f,
                const GeometricField<Type,fvsPatchField,surfaceMesh>& fSub
            ) const;

            //- Copy a surface field defined on the mesh to the subMesh
            template<class Type>
            void copyMeshToSubMesh
            (
                const GeometricField<Type,fvsPatchField,surfaceMesh>& f,
                GeometricField<Type,fvsPatchField,surfaceMesh>& fSub
            ) const;

        // Write

            //- Write class output
            void write(Ostream&) const;


};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshToSubMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
