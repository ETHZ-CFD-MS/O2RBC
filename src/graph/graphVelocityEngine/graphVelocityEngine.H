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
    Foam::graphVelocityEngine

Description
    Compute divergence-free advective fluxes through a tube network. The fluxes
    are computed from velocities defined on graph edges.

    The velocities are interpolated onto the lumen mesh and then projected to
    form a set of divergence-free fluxes.

    This class supports lumen meshes with multiple connected components

SourceFiles
    graphVelocityEngine.C

\*---------------------------------------------------------------------------*/

#ifndef graphVelocityEngine_H
#define graphVelocityEngine_H

#include "fvCFD.H"
#include "geometricEdgeGraph.H"
#include "vascularGraphRegions.H"

#include "vectorList.H"
#include "meshToSubMesh.H"
#include "piecewiseConstantInterpolationTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;

// Forward declaration of friend functions and operators
class graphVelocityEngine;
Ostream& operator<<(Ostream&, const graphVelocityEngine&);


/*---------------------------------------------------------------------------*\
                         Class graphVelocityEngine Declaration
\*---------------------------------------------------------------------------*/

class graphVelocityEngine
:
    public IOdictionary
{
    // Private data

        //- Fixed Eulerian mesh
        const fvMesh& mesh_;

        //- Submesh for lumen
        fvMesh lumenMesh_;

        //- Graph structure
        const geometricEdgeGraph& graph_;

        //- Vascular graph regions
        const vascularGraphRegions& vesselRegions_;

        //- Flux field defined on the whole mesh
        surfaceScalarField phi_;

        //- Flux field defined on the lumen submesh
        surfaceScalarField phiLumen_;

        //- Scalar field for the flux correction
        volScalarField fCorr_;

        //- Velocity field mapped from graph velocities
        volVectorField UMappedLumen_;

        //- Velocity field reconstructed from phi_
        volVectorField UReconstructed_;

        //- Time-dependent velocity data on each edge
        //  Uses a single-linked list so that binary data can be read
        piecewiseConstantInterpolationTable<SLList<scalar> > edgeVelocities_;

        //- Name of the patch that cuts the lumen
        const word lumenPatchName_;

        //- Name of the patch for the capillary wall
        const word wallPatchName_;

        //- Object for the copy of fields between mesh_ and lumenMesh_.
        meshToSubMesh meshToLumen_;

        //- Reference cell for the solution to fCorr.
        label fCorrRefCell_;

        //- Reference value for the solution to fCorr.
        scalar fCorrRefValue_;


    // Private Member Functions

        //- From the velocity defined on the graph edges, compute advective
        //  fluxes on the lumen mesh.
        void transformGraphVelocityToLumenFluxes();

        //- From the velocity defined on the graph edges, compute the
        // advection velocity on the lumen mesh.
        void transformGraphVelocityToLumenVelocity();

        //- Compute the advection flux inside the domain, given edge velocities and
        //  face labels in the lumen mesh.
        void computeInteriorFlux
        (
            const vectorList& edgeVelocityVectors,
            const labelField& faceEdgeLabelsLumen
        );

        //- Compute the advection flux on a patch, given edge velocities and
        //  face labels in the lumen mesh.
        void computePatchFlux
        (
            const fvPatch&,
            const vectorList& edgeVelocityVectors,
            const labelField& faceEdgeLabelsLumen
        );

        //- Compute the fluxes at processor boundaries
        void computeProcessorFluxes();

        //- Correct the boundary fluxes so that the inflow equals the outflow.
        void correctNetBoundaryFlux();

        //- Correct the boundary fluxes so that the inflow equals the outflow on
        //  each cellZone. This can be used when the mesh has multiple connected
        //  components with corresponding cell zones.
        void correctCellZoneNetBoundaryFlux();

        //- Synchronize value of field across processor boundaries
        //- TODO: move this somewhere else, or make template out of it.
        // void synchronizeProcessorPatchField(GeometricField<Type, PatchField, GeoMesh>&);
        void synchronizeProcessorPatchField(volVectorField&);

        //- Return the patch ID for the capillary wall
        label wallPatchID() const;

        //- Return the patch ID for the lumen
        label lumenPatchID() const;

        //- Return the patch for the capillary wall
        const fvPatch& wallPatch() const;

        //- Return the patch that cuts the lumen
        const fvPatch& lumenPatch() const;

        //- Return the zoneID of the cell adjacent to face patch
        inline label patchFaceToCellZoneID
        (
            const label patchFaceI, 
            const fvPatch&
        ) const;

        //- Disallow default bitwise copy construct
        graphVelocityEngine(const graphVelocityEngine&);

        //- Disallow default bitwise assignment
        void operator=(const graphVelocityEngine&);


public:

    //- Runtime type information
    TypeName("graphVelocityEngine");


    // Constructors

        //- Construct from components
        graphVelocityEngine
        (
            const fvMesh& mesh,
            const word& lumenName,
            const geometricEdgeGraph&,
            const vascularGraphRegions&
        );


    //- Destructor
    ~graphVelocityEngine();


    // Member Functions

        const surfaceScalarField& conservativePhi();

        // Write

            //- WriteData member function required by regIOobject.
            bool writeData(Ostream& os) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const graphVelocityEngine&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "graphVelocityEngineI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
