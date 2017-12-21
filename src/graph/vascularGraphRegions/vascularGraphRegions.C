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

#include "vascularGraphRegions.H"

#include "circularTubeGraph.H"
#include "meshToMeshMoving.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vascularGraphRegions, 1);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool
Foam::vascularGraphRegions::edgeOverlapsMesh(const label edgeIndex) const
{
    pointField tubeSkeleton = graph_.tube(edgeIndex).skeletonPoints();
    const boundBox& meshBb = mesh_.bounds();

    forAll(tubeSkeleton, i)
    {
        if (meshBb.contains(tubeSkeleton[i]))
        {
            return true;
        }
    }

    return false;
}


const Foam::labelList
Foam::vascularGraphRegions::overlappingEdgeIndices() const
{
    labelList indexList;

    forAll(graph_.edgeIndices(), eI)
    {
        if (edgeOverlapsMesh(eI))
        {
            indexList.append(eI);
        }
    }

    return indexList;
}


void
Foam::vascularGraphRegions::interpolateRegions
(
    const volScalarField& inLumenCylinder,
    const volScalarField& inWallCylinder,
    const volScalarField& normalizedAxisCoord
 )
{
    inLumen_              = 0.0;
    inWall_               = 0.0;
    inTissue_             = 0.0;
    cellEdgeLabels_       = labelField(mesh_.nCells(), -1);
    faceEdgeLabels_       = labelField(mesh_.nFaces(), -1);
    edgeCurvilinearCoord_ = 0.0;

    volScalarField oneOnCylinder(inLumenCylinder);
    oneOnCylinder   = 1.0;

    volScalarField edgeIndicator(inLumen_);
    volScalarField edgeLocalInLumen(inLumen_);
    volScalarField cylinderAxisNormalizedCoord(inLumen_);

    const labelList overlappingEdges = overlappingEdgeIndices();

    forAll(tubeMeshes_, i)
    {
        label edgeI = overlappingEdges[i];

        // create interpolation object
        meshToMeshMoving interp
        (
            tubeMeshes_[i],
            mesh_,
            Foam::meshToMeshMoving::imCellRelativeVolumeWeight,
            /* interpAllPatches */ true,
            /* serialInterp */ true
        );

        // Interpolation
        // The boolean argument indicates not to use a weighting factor for the old target values.
        // This allows the results for vascular graph regions to be correct in
        // Eulerian grid cells where cylinders overlap.
        interp.mapSrcToTgt
        (
            inLumenCylinder.internalField(), 
            plusEqOp<double>(), 
            inLumen_.primitiveFieldRef(), 
            false 
        );

        interp.mapSrcToTgt
        (
            inWallCylinder.internalField(), 
            plusEqOp<double>(), 
            inWall_.primitiveFieldRef(),
            false
        );

        // Compute cell labels and edge curvilinear coordinates
        cylinderAxisNormalizedCoord = 0.0;
        interp.mapSrcToTgt
        (
            normalizedAxisCoord.internalField(), 
            plusEqOp<double>(), 
            cylinderAxisNormalizedCoord.primitiveFieldRef(),
            false
        );

        edgeIndicator = 0.0;
        interp.mapSrcToTgt
        (
            oneOnCylinder.internalField(), 
            plusEqOp<double>(), 
            edgeIndicator.primitiveFieldRef(), 
            false 
        );

        forAll(mesh_.cells(), i)
        {
            if (edgeIndicator[i] > 0)
            {
                scalar edgeLength = graph_.edge(edgeI).length();
                cellEdgeLabels_[i] = edgeI;
                // divide the interpolated normalized coordinate by the value of
                // edgeIndicator to get the actual normalized coordinate.
                edgeCurvilinearCoord_[i] =   edgeLength
                                           * cylinderAxisNormalizedCoord[i]
                                           / edgeIndicator[i];
            }
        }
    }

    // Compute the face edge labels by using the edgeCellLabels
    forAll(mesh_.faces(), faceI)
    {
        label ownerEdgeLabel   = cellEdgeLabels_[mesh_.faceOwner()[faceI]];
        faceEdgeLabels_[faceI] = ownerEdgeLabel;
    }

    // Convert the integer edge labels to scalar edge labels
    forAll(mesh_.cells(), i)
    {
        cellEdgeScalarLabels_.primitiveFieldRef()[i] = cellEdgeLabels_[i];
    }
    for(int faceI = 0; faceI < mesh_.nInternalFaces() ; faceI++)
    {
        faceEdgeScalarLabels_.primitiveFieldRef()[faceI] = faceEdgeLabels_[faceI];
    }

    // Ensure that the sum of inWall and inLumen does not exceed one.
    // If it is the case, clamp inLumen and inWall to a maximal value of one.
    // Then, further clamp inWall so that inLumen + inWall <= 1.0, to remove 
    // disconnected wall region parts that may arise at bifurcations.
    forAll(mesh_.cells(), i)
    {
        inLumen_[i] = min(inLumen_[i], 1.0);
        inWall_ [i] = min(inWall_ [i], 1.0 - inLumen_[i]);
    }

    // compute inTissue
    inTissue_ = 1.0 - inLumen_ - inWall_;

    // set values on patch to be equal to the neighboring internal field.
    forAll(mesh_.boundary(), i)
    {
        inLumen_ .boundaryFieldRef()[i] = inLumen_ .boundaryField()[i].patchInternalField();
        inWall_  .boundaryFieldRef()[i] = inWall_  .boundaryField()[i].patchInternalField();
        inTissue_.boundaryFieldRef()[i] = inTissue_.boundaryField()[i].patchInternalField();

        cellEdgeScalarLabels_.boundaryFieldRef()[i] 
                = cellEdgeScalarLabels_.boundaryField()[i].patchInternalField();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vascularGraphRegions::vascularGraphRegions
(
    const fvMesh& mesh, 
    const circularTubeGraph& graph
)
:
    mesh_(mesh),
    graph_(graph),
    tubeMeshes_(graph_.nEdges()),
    inLumen_
    (
        IOobject
        (
            "inLumen",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    inWall_
    (
        IOobject
        (
            "inWall",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    inTissue_
    (
        IOobject
        (
            "inTissue",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    edgeCurvilinearCoord_
    (
        IOobject
        (
            "edgeCurvilinearCoord",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    cellEdgeLabels_
    (
        IOobject
        (
            "cellEdgeLabels",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_.nCells()
    ),
    cellEdgeScalarLabels_
    (
        IOobject
        (
            "cellEdgeScalarLabels",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    ),
    faceEdgeLabels_
    (
        IOobject
        (
            "faceEdgeLabels",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_.nFaces()
    ),
    faceEdgeScalarLabels_
    (
        IOobject
        (
            "faceEdgeScalarLabels",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    )
{
    // convert the scalar labels to integer labels
    forAll(mesh.cells(), i)
    {
        cellEdgeLabels_[i] = (int) round(cellEdgeScalarLabels_.internalField()[i]);
    }

    // read the edge meshes
    const labelList overlappingEdges = overlappingEdgeIndices();

    forAll(overlappingEdges, i)
    {
        label edgeI = overlappingEdges[i];

        tubeMeshes_.append
        (
            new polyMesh
            (
                IOobject
                (
                    tubeMapperVisitor::tubeMeshName(edgeI),
                    mesh_.time().constant(),
                    mesh_.time(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false // do not register in database
                ),
                xferCopy(pointField()),
                xferCopy(faceList()),
                xferCopy(labelList()),
                xferCopy(labelList()),
                false // no parallel comms
            )
        );
    }
}


Foam::vascularGraphRegions::vascularGraphRegions
(
    const fvMesh& mesh, 
    const circularTubeGraph& graph,
    const word& cylinderMeshName
)
:
    mesh_(mesh),
    graph_(graph),
    tubeMeshes_(),
    inLumen_
    (
        IOobject
        (
            "inLumen",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimless
    ),
    inWall_
    (
        IOobject
        (
            "inWall",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimless
    ),
    inTissue_
    (
        IOobject
        (
            "inTissue",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimless
    ),
    edgeCurvilinearCoord_
    (
        IOobject
        (
            "edgeCurvilinearCoord",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimless
    ),
    cellEdgeLabels_
    (
        IOobject
        (
            "cellEdgeLabels",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_.nCells()
    ),
    cellEdgeScalarLabels_
    (
        IOobject
        (
            "cellEdgeScalarLabels",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimless
    ),
    faceEdgeLabels_
    (
        IOobject
        (
            "faceEdgeLabels",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_.nFaces()
    ),
    faceEdgeScalarLabels_
    (
        IOobject
        (
            "faceEdgeScalarLabels",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimless
    )
{
    calculate(cylinderMeshName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vascularGraphRegions::~vascularGraphRegions()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
Foam::vascularGraphRegions::calculate(const word& cylinderMeshName)
{
    Info<< "Constructing cylindrical meshes for graph tubes..." << endl;

    tubeMapperVisitor tubeMapper
    (
        IOobject
        (
            cylinderMeshName,
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::MUST_READ
        )
    );

    // Construct graph and read required fields
    volScalarField inLumenCylinder
    (
        IOobject
        (
            "inLumenCylinder",
            tubeMapper.cylinderMesh().time().caseConstant(),
            tubeMapper.cylinderMesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        tubeMapper.cylinderMesh()
    );
    volScalarField inWallCylinder
    (
        IOobject
        (
            "inWallCylinder",
            tubeMapper.cylinderMesh().time().caseConstant(),
            tubeMapper.cylinderMesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        tubeMapper.cylinderMesh()
    );

    const labelList overlappingEdges = overlappingEdgeIndices();

    forAll(overlappingEdges, i)
    {
        label edgeI = overlappingEdges[i];

        if (debug)
        {
            Pout<< "Foam::vascularGraphRegions::calculate(const word&): "
                << "Constructing mesh for edge " << edgeI << endl;
        }

        // Pass the tube mapper to current tube.
        // This will call tubeMapper.visit(tube) and compute the point
        // positions. See the Visitor pattern.
        graph_.tube(edgeI).accept(tubeMapper);

        tubeMeshes_.append
        (
            tubeMapper.mappedMesh(edgeI)
        );
    }

    autoPtr<volScalarField> pNormalizedAxisCoord 
                    = tubeMapper.normalizedCylinderAxisCoord();
    const volScalarField& normalizedAxisCoord = pNormalizedAxisCoord();

    interpolateRegions(inLumenCylinder, inWallCylinder, normalizedAxisCoord);
}


void
Foam::vascularGraphRegions::writeFields() const
{
    inLumen_.write();
    inWall_.write();
    inTissue_.write();
    cellEdgeLabels_.write();
    cellEdgeScalarLabels_.write();
    faceEdgeLabels_.write();
    faceEdgeScalarLabels_.write();
    edgeCurvilinearCoord_.write();
}


void
Foam::vascularGraphRegions::writeMeshes() const
{
    forAll(tubeMeshes_, i)
    {
        tubeMeshes_[i].write();
    }
}



// ************************************************************************* //
