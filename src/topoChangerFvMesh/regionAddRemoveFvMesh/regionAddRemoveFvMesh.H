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

Class
    Foam::regionAddRemoveFvMesh

Description
    Add and remove regions from an original mesh to create a single new mesh.

    Large parts of the code were taken from the class mergePolyMesh.

SourceFiles
    regionAddRemoveFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef regionAddRemoveFvMesh_H
#define regionAddRemoveFvMesh_H

#include "topoChangerFvMesh.H"
#include "polyTopoChange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class regionAddRemoveFvMesh Declaration
\*---------------------------------------------------------------------------*/

class regionAddRemoveFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Topological change to accumulated all mesh changes
        mutable autoPtr<polyTopoChange> meshModPtr_;

        //- Patch names
        DynamicList<word> patchNames_;

        //- Patch dictionaries
        DynamicList<dictionary> patchDicts_;

        //- Point zone names
        DynamicList<word> pointZoneNames_;

        //- Face zone names
        DynamicList<word> faceZoneNames_;

        //- Cell zone names
        DynamicList<word> cellZoneNames_;

        //- Map for fields after topological changes
        autoPtr<mapPolyMesh> meshMap_;


    // Private Member Functions

        // Return the polyTopoChange object
        polyTopoChange& meshModifier() const;

        //- Return patch index given a name and type
        label patchIndex(const polyPatch&);

        //- Return zone index given a list of active zones and a name
        label zoneIndex(DynamicList<word>&, const word&);

        //- Change the patch ID of a face
        void changePatchID
        (
            const label faceID,
            const label patchID
        );

        //- Disallow default bitwise copy construct
        regionAddRemoveFvMesh(const regionAddRemoveFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const regionAddRemoveFvMesh&);


public:

    //- Runtime type information
    TypeName("regionAddRemoveFvMesh");


    // Constructors

        //- Construct from IOobject
        regionAddRemoveFvMesh(const IOobject&);


    //- Destructor
    ~regionAddRemoveFvMesh();


    // Member Functions

        //- Access the mapPolyMesh object produced by topological changes
        autoPtr<mapPolyMesh> meshMap()
        {
            return meshMap_;
        }

        //- Add a mesh
        void addMesh
        (
            const polyMesh& m
        );

        //- Remove a zone
        void removeMesh
        (
            const word& cellZoneName,
            const word& pointZoneName,
            const word& faceZoneName
        );

        //- Merge meshes
        void merge();

        //- Filter out the empty patches
        void filterPatches();
        
        //- Update the mesh for topological changes
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

