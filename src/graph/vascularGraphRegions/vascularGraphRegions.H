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
    Foam::vascularGraphRegions

Description
    Computes the vascular regions (plasma, capillary wall and tissue) for an
    instance of circularTubeGraph. 
    
    This computation is performed by constructing a list of tube meshes (one 
    per capillary) and interpolating the indicator fields (inLumenCylinder and 
    inWallCylinder) onto the mesh used for oxygen transport.

SourceFiles
    vascularGraphRegions.C

\*---------------------------------------------------------------------------*/

#ifndef vascularGraphRegions_H
#define vascularGraphRegions_H

#include "fvCFD.H"

#include "tubeMapperVisitor.H"
#include "labelIOField.H"

namespace Foam
{

// Forward declaration of classes
class Ostream;
class circularTubeGraph;

// Forward declaration of friend functions and operators
class vascularGraphRegions;
Ostream& operator<<(Ostream&, const vascularGraphRegions&);


/*---------------------------------------------------------------------------*\
                         Class vascularGraphRegions Declaration
\*---------------------------------------------------------------------------*/

class vascularGraphRegions
{
    // Private data

        //- Reference to Eulerian mesh
        const fvMesh& mesh_;

        //- Reference to graph
        const circularTubeGraph& graph_;

        //- Collection of meshes for the graph tubes
        PtrList<polyMesh> tubeMeshes_;

        //- Field for lumen volume fraction on mesh_
        volScalarField inLumen_;

        //- Field for capillary wall volume fraction on mesh_
        volScalarField inWall_;

        //- Field for tissue volume fraction on mesh_
        volScalarField inTissue_;

        //- Field for the edge curvilinear coordinate
        volScalarField edgeCurvilinearCoord_;

        //- Cell-centered labels of the graph edges.
        //  The internal edge indices are used (between 0 and IOgraph::nEdges())
        labelIOField cellEdgeLabels_;
        
        //- Cell-centered labels of the graph edges (scalar valued)
        //  Use for input/output since decomposePar does not decompose label
        //  fields.
        volScalarField cellEdgeScalarLabels_;

        //- Face-centered labels of the graph edges.
        //  The internal edge indices are used (between 0 and IOgraph::nEdges())
        labelIOField faceEdgeLabels_;

        //- Face-centered labels of the graph edges (scalar valued)
        //  Use for input/output since decomposePar does not decompose label
        //  fields.
        surfaceScalarField faceEdgeScalarLabels_;

    // Private Member Functions
    
        //- Return whether the given edge overlaps the Eulerian mesh.
        //  Takes a local index as argument.
        bool edgeOverlapsMesh(const label edgeIndex) const;

        //- List of edges that overlap the Eulerian mesh.
        //  Returns local indices.
        const labelList overlappingEdgeIndices() const;

        //- Interpolate the region fields and the normalized cylinder axis
        //  coordinate from the base cylinder mesh to the Eulerian mesh.
        void interpolateRegions
        (
            const volScalarField& inLumenCylinder,
            const volScalarField& inWallCylinder,
            const volScalarField& normalizedAxisCoord
        );

        //- Disallow default bitwise copy construct
        vascularGraphRegions(const vascularGraphRegions&);

        //- Disallow default bitwise assignment
        void operator=(const vascularGraphRegions&);


public:

    // Runtime type information
    ClassName("vascularGraphRegions");

    // Constructors

        //- Construct by reading the vascular regions and the associated edge
        //  meshes
        vascularGraphRegions
        (
            const fvMesh& mesh, 
            const circularTubeGraph& 
        );

        //- Construct by calculating the vascular regions
        vascularGraphRegions
        (
            const fvMesh& mesh, 
            const circularTubeGraph&, 
            const word& cylinderMeshName
        );


    //- Destructor
    ~vascularGraphRegions();


    // Member Functions

        //- Calculate regions
        void calculate(const word& cylinderMeshName);

        //- Constant access to volume fraction of capillary lumen
        const volScalarField& inLumen() const
        {
            return inLumen_;
        }
        
        //- Constant access to volume fraction of capillary wall
        const volScalarField& inWall() const
        {
            return inWall_;
        }

        //- Constant access to volume fraction of tissue
        const volScalarField& inTissue() const
        {
            return inTissue_;
        }

        //- Constant access to cell edge labels
        const labelField& cellEdgeLabels() const
        {
            return cellEdgeLabels_;
        }

        //- Constant access to face edge labels
        const labelField& faceEdgeLabels() const
        {
            return faceEdgeLabels_;
        }

        //- Constant access to the edge curvilinear coordinate
        const volScalarField& edgeCurvilinearCoord() const
        {
            return edgeCurvilinearCoord_;
        }

        //- Write the fields computed by this class
        void writeFields() const;

        //- Write the meshes computed by this class
        void writeMeshes() const;

    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const vascularGraphRegions&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
