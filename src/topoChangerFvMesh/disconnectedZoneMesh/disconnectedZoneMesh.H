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
    Foam::disconnectedZoneMesh

Description
    Class for a mesh that consists of disconnected zones. Zone motion is 
    supported. Setting and getting of fields on individual zones is 
    implemented.
    
    It is assumed that the cell, point and face zones are 
    consistent, that is, these zones correspond to a given connected component 
    of the mesh.

SourceFiles
    disconnectedZoneMesh.C
    disconnectedZoneMeshTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef disconnectedZoneMesh_H
#define disconnectedZoneMesh_H

#include "regionAddRemoveFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class disconnectedZoneMesh Declaration
\*---------------------------------------------------------------------------*/

class disconnectedZoneMesh
:
    public regionAddRemoveFvMesh
{
    // Private data


    // Private Member Functions


        //- Disallow default bitwise copy construct
        disconnectedZoneMesh(const disconnectedZoneMesh&);

        //- Disallow default bitwise assignment
        void operator=(const disconnectedZoneMesh&);


public:

    //- Runtime type information
    TypeName("disconnectedZoneMesh");


    // Constructors

        //- Construct from IOobject
        disconnectedZoneMesh(const IOobject&);


    //- Destructor
    ~disconnectedZoneMesh();


    // Member Functions

        // Access
        
            //- Number of cells in a given zone
            label nZoneCells(const label zoneI) const;

            //- List of patch indices and boundary face indices of a given zone.
            //  The first element of the tuple is the patch index, the second is
            //  the face index (local to the patch).
            labelPairList zoneBoundaryFaces(const label zoneI) const;

            //- Constant access to field values on a given zone
            template<class Type>
            Field<Type> getFieldOnZone
            (
                const Field<Type>&,
                const label zoneI
            ) const;

            //- Constant access to geometric field values on a given zone
            template<class Type>
            Field<Type> getFieldOnZone
            (
                const GeometricField<Type,fvPatchField,volMesh>&, 
                const label zoneI
            ) const;

            //- Bounding box of a zone
            boundBox zoneBounds(const label zoneI) const;


        // Check

            //- Whether a zone index has associated non-empty zones
            bool hasNonEmptyZone(const label zoneI) const;


        // Edit

            //- Set field values on a given zone and its boundary. The argument
            //  zoneField should have the size of the zone.
            template<class Type>
            void setFieldOnZone
            (
                Field<Type>& f, 
                const Field<Type>& zoneField,
                const label zoneI
            ) const;

            //- Set field to a constant value on a given zone
            template<class Type>
            void setFieldOnZone
            (
                Field<Type>& f, 
                const Type value,
                const label zoneI
            ) const;

            //- Set geometric field to a constant value on a given zone and its boundary
            template<class Type>
            void setFieldOnZone
            (
                GeometricField<Type,fvPatchField,volMesh>& f, 
                const Type value,
                const label zoneI
            ) const;

            //- Set geometric field to a constant value on a zone boundary
            template<class Type>
            void setFieldOnZoneBoundary
            (
                GeometricField<Type,fvPatchField,volMesh>& f, 
                const Type value,
                const label zoneI
            ) const;

            //- Set boundary field to the patch internal field on a given zone
            //  boundary
            template<class Type>
            void setBoundaryFieldToPatchInternalFieldOnZone
            (
                GeometricField<Type,fvPatchField,volMesh>& f, 
                const label zoneI
            ) const;

            //- Move the points of a given zone
            void movePoints(const pointField&, const label zoneI);

        // Write


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
    #include "disconnectedZoneMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
