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

#include "graphVelocityEngine.H"

#include "processorFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(graphVelocityEngine, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
Foam::graphVelocityEngine::transformGraphVelocityToLumenVelocity()
{
    const labelField& cellEdgeLabels = vesselRegions_.cellEdgeLabels();
    const volScalarField& edgeCoord = vesselRegions_.edgeCurvilinearCoord();

    const scalar time = mesh_.time().timeOutputValue();
    const scalarList& currentEdgeVelocities = scalarList(edgeVelocities_(time));

    scalar maxEdgeAbsVelocity = 0.0;

    labelField cellEdgeLabelsLumen(lumenMesh_.nCells(), -1);
    scalarField edgeCoordLumen(lumenMesh_.nCells(), -1.0);

    meshToLumen_.copyMeshToSubMeshCells(cellEdgeLabels, cellEdgeLabelsLumen);
    meshToLumen_.copyMeshToSubMeshCells(edgeCoord.internalField(), edgeCoordLumen);

    // Internal field
    forAll(UMappedLumen_, cellI)
    {
        label edgeI = cellEdgeLabelsLumen[cellI];
        if (edgeI >= 0)
        {
            const scalar edgeVelocity = currentEdgeVelocities[edgeI];
            const scalar sCoord = edgeCoordLumen[cellI];
            const vector flowDirection = graph_.edge(edgeI).tangentVector(sCoord);

            UMappedLumen_[cellI] = edgeVelocity*flowDirection;
            maxEdgeAbsVelocity = std::max(edgeVelocity, maxEdgeAbsVelocity);
        }
    }

    if (debug)
    {
        Info<< "Foam::graphVelocityEngine::transformGraphVelocityToLumenVelocity:: "
            << "Maximal velocity along an edge: " << maxEdgeAbsVelocity << endl;
    }

    // Patch fields
    const polyBoundaryMesh& patches = lumenMesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        UMappedLumen_.boundaryFieldRef()[patchI] = 
            UMappedLumen_.boundaryField()[patchI].patchInternalField();
    }

    UMappedLumen_.boundaryFieldRef()[wallPatchID()] = vector::zero;
}


void
Foam::graphVelocityEngine::computeProcessorFluxes()
{
    const fvBoundaryMesh& patches = lumenMesh_.boundary();

    forAll(patches, patchI)
    {
        const fvPatch& patch = patches[patchI];
        const vectorField& patchNormalVectors = patch.Sf();

        fvsPatchField<scalar>& phiPf = phiLumen_.boundaryFieldRef()[patchI];
        const fvPatchField<vector>& UPf = UMappedLumen_.boundaryField()[patchI];

        if (isA<processorFvPatch>(patch))
        {
            forAll(patch, i)
            {
                phiPf[i] = UPf[i] & patchNormalVectors[i];
            }
        }
    }
}


void
Foam::graphVelocityEngine::correctNetBoundaryFlux()
{
    // sum up the fluxes on the boundary
    // TODO: if needed, introduced a more refined version that first adapts the
    // fluxes for faces that touch the wall (i.e., cells with inWall_ > 0).
    const label patchID = lumenPatchID();
    const fvPatch& cPatch = lumenPatch();

    scalar fluxSum    = gSum(phiLumen_.boundaryField()[patchID]);
    scalar surfaceSum = gSum(cPatch.magSf());

    if (debug)
    {
        Info<< "Foam::graphVelocityEngine::correctNetBoundaryFlux:: "
            << "Sum of the fluxes over lumen patch = " << fluxSum << ", "
            << "Sum of surface element area = " << surfaceSum << endl;
    }

    // correct the fluxes
    forAll(cPatch, i)
    {
        // compute a surface-weighted fraction of fluxSum
        scalar fluxCorr = fluxSum*cPatch.magSf()[i]/surfaceSum;
        phiLumen_.boundaryFieldRef()[patchID][i] -= fluxCorr;
    }

    if (debug)
    {
        fluxSum = gSum(phiLumen_.boundaryField()[patchID]);
        Info<< "Foam::graphVelocityEngine::correctNetBoundaryFlux:: "
            << "Corrected flux sum = " << fluxSum << endl;
    }
}


void
Foam::graphVelocityEngine::correctCellZoneNetBoundaryFlux()
{
    const label patchID = lumenPatchID();
    const fvPatch& lumenInOutPatch = lumenPatch();
    const label nZones = lumenMesh_.cellZones().size();

    // sum up the fluxes on the boundary
    scalarList fluxSum(nZones, 0.0);
    scalarList surfaceSum(nZones, 0.0);

    forAll(lumenInOutPatch, i)
    {
        const label zoneI = patchFaceToCellZoneID(i, lumenInOutPatch);

        fluxSum[zoneI] += phiLumen_.boundaryField()[patchID][i];
        surfaceSum[zoneI] += lumenInOutPatch.magSf()[i];
    }

    reduce(fluxSum,    sumOp<scalarList>());
    reduce(surfaceSum, sumOp<scalarList>());

    if (debug)
    {
        Info<< "Foam::graphVelocityEngine::correctCellZoneNetBoundaryFlux:: "
            << "Sum of the fluxes over lumen patch = " << fluxSum << ", "
            << "Sum of surface element area = " << surfaceSum << endl;
    }

    // correct the fluxes
    forAll(lumenInOutPatch, i)
    {
        const label zoneI = patchFaceToCellZoneID(i, lumenInOutPatch);

        // compute a surface-weighted fraction of the flux sum
        scalar fluxCorr = fluxSum[zoneI]
                        * lumenInOutPatch.magSf()[i]/surfaceSum[zoneI];

        phiLumen_.boundaryFieldRef()[patchID][i] -= fluxCorr;
    }

    // check the results
    if (debug)
    {
        scalarList correctedFluxSum(nZones, 0.0);

        forAll(lumenInOutPatch, i)
        {
            const label zoneI = patchFaceToCellZoneID(i, lumenInOutPatch);

            correctedFluxSum[zoneI] += phiLumen_.boundaryField()[patchID][i];
        }

        reduce(correctedFluxSum, sumOp<scalarList>());

        Info<< "Foam::graphVelocityEngine::correctCellZoneNetBoundaryFlux:: "
            << "Corrected flux sum = " << correctedFluxSum << endl;
    }
}


void
Foam::graphVelocityEngine::synchronizeProcessorPatchField
(
    volVectorField& vf
)
{
    const polyBoundaryMesh& patches = vf.mesh().boundaryMesh();

    FieldField<fvPatchField, vector> oldBoundaryField = vf.boundaryField();

    forAll(patches, patchI)
    {
        vf.boundaryFieldRef()[patchI].initEvaluate(Pstream::blocking);
    }
    
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        fvPatchField<vector>& pf = vf.boundaryFieldRef()[patchI];
        // The following call puts the values of the patch neighbour field
        // info pf.
        pf.evaluate(Pstream::blocking);

        if (isA<processorPolyPatch>(pp))
        {
            // TODO: use delta coefficients on each side instead
            pf = 0.5*(pf.patchNeighbourField() + oldBoundaryField[patchI]);
        }
    }
}


Foam::label
Foam::graphVelocityEngine::wallPatchID() const
{
    return lumenMesh_.boundary().findPatchID(wallPatchName_);
}


Foam::label
Foam::graphVelocityEngine::lumenPatchID() const
{
    return lumenMesh_.boundary().findPatchID(lumenPatchName_);
}


const Foam::fvPatch& 
Foam::graphVelocityEngine::wallPatch() const
{
    return lumenMesh_.boundary()[wallPatchID()];
}


const Foam::fvPatch& 
Foam::graphVelocityEngine::lumenPatch() const
{
    return lumenMesh_.boundary()[lumenPatchID()];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::graphVelocityEngine::graphVelocityEngine
(
    const fvMesh& mesh,
    const word& lumenName,
    const geometricEdgeGraph& graph,
    const vascularGraphRegions& vesselRegions
)
:
    IOdictionary
    (
        IOobject
        (
            "graphVelocityEngineDict",
            mesh.time().caseSystem(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    lumenMesh_
    (
        IOobject
        (
            lumenName,
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    ),
    graph_(graph),
    vesselRegions_(vesselRegions),
    phi_
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("phi", dimVolume/dimTime, 0.0)
    ),
    phiLumen_
    (
        IOobject
        (
            "phiLumen",
            lumenMesh_.time().timeName(),
            lumenMesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        lumenMesh_,
        dimensionedScalar("phiLumen", dimVolume/dimTime, 0.0)
    ),
    fCorr_
    (
        IOobject
        (
            "fCorr",
            lumenMesh_.time().timeName(),
            lumenMesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        lumenMesh_
    ),
    UMappedLumen_
    (
        IOobject
        (
            "UMappedLumen",
            lumenMesh_.time().timeName(),
            lumenMesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        lumenMesh_,
        dimensionedVector("UMappedLumen", dimVelocity, vector::zero)
    ),
    UReconstructed_
    (
        IOobject
        (
            "UReconstructed",
            lumenMesh_.time().timeName(),
            lumenMesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        lumenMesh_,
        dimensionedVector("UReconstructed", dimVelocity, vector::zero)
    ),
    edgeVelocities_
    (
        IOdictionary
        (
            IOobject
            (
                "edgeVelocities",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            )
        ),
        "edgeVelocities"
    ),
    lumenPatchName_(lookup("lumenPatch")),
    wallPatchName_(lookup("wallPatch")),
    meshToLumen_(mesh_, lumenMesh_),
    fCorrRefCell_(0),
    fCorrRefValue_(0.0)
{
    // sanity checks
    forAll(edgeVelocities_, i)
    {
        if (edgeVelocities_[i].second().size() != graph_.nEdges())
        {
            FatalErrorIn
            (
                "graphVelocityEngine(const fvMesh&, const word&, "
                "const geometricEdgeGraph&, const vascularGraphRegions&)"
            ) << "The number of elements in edgeVelocities[" << i << "] "
              << "does not equal the number of edges in the graph ("
              << graph_.nEdges() << ")." << nl
              << abort(FatalError);
        }
    }

    if (lumenPatchID() == -1)
    {
        FatalErrorIn
        (
            "graphVelocityEngine(const fvMesh&, const word&, "
            "const geometricEdgeGraph&, const vascularGraphRegions&)"
        ) << "Patch " << lumenPatchName_ << " not found."
          << abort(FatalError);
    }

    if (wallPatchID() == -1)
    {
        FatalErrorIn
        (
            "graphVelocityEngine(const fvMesh&, const word&, "
            "const geometricEdgeGraph&, const vascularGraphRegions&)"
        ) << "Patch " << wallPatchName_ << " not found."
          << abort(FatalError);
    }

    if (lumenMesh_.cellZones().size() == 0)
    {
        WarningIn
        (
            "graphVelocityEngine(const fvMesh&, const word&, "
            "const geometricEdgeGraph&, const vascularGraphRegions&)"
        ) << "The lumen mesh contains no cell zones. Therefore, this class "
          << "will assume that the lumen mesh consists of only one connected "
          << "component." << endl;

        autoPtr<cellZone> czPtr
        (
            new cellZone
            (
                "region0", 
                identity(lumenMesh_.nCells()), 
                0, 
                lumenMesh_.cellZones()
            )
        );
        List<cellZone*>  czList(1, czPtr.ptr());
        List<pointZone*> pzList;
        List<faceZone*>  fzList;
        lumenMesh_.addZones(pzList, fzList, czList);
    }

    // read reference value and cell
    setRefCell
    (
        fCorr_, 
        lumenMesh_.solutionDict().subDict("graphVelocityEngine"), 
        fCorrRefCell_, 
        fCorrRefValue_
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::graphVelocityEngine::~graphVelocityEngine()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::surfaceScalarField&
Foam::graphVelocityEngine::conservativePhi()
{
    scalar currentTime  = mesh_.time().timeOutputValue();
    scalar previousTime = currentTime - mesh_.time().deltaTValue();

    // if velocity changed, or if it is the first time step
    if (edgeVelocities_.valueChanged(previousTime, currentTime)
     || mesh_.time().timeIndex() <= mesh_.time().startTimeIndex() + 1)
    {
        transformGraphVelocityToLumenVelocity();
        synchronizeProcessorPatchField(UMappedLumen_);

        phiLumen_ = fvc::interpolate(UMappedLumen_) & lumenMesh_.Sf();
        computeProcessorFluxes();
        correctCellZoneNetBoundaryFlux();
        phiLumen_.boundaryFieldRef()[wallPatchID()] = 0.0;

        if (debug)
        {
            Info<< "Foam::graphVelocityEngine::conservativePhi(): "
                << "Maximum of divergence of phiLumen before projection: "
                << max(mag(fvc::div(phiLumen_))) << endl;
        }

        // make phiLumen divergence-free
        fvScalarMatrix fCorrEqn
        (
            fvm::laplacian(fCorr_) == fvc::div(phiLumen_)
        );

        fCorrEqn.setReference(fCorrRefCell_, fCorrRefValue_);
        fCorrEqn.solve();

        phiLumen_ -= fCorrEqn.flux();

        // copy phiLumen from the lumen mesh to the mesh
        meshToLumen_.copySubMeshToMesh(phi_, phiLumen_);

        if (debug)
        {
            Info<< "Foam::graphVelocityEngine::conservativePhi(): "
                << "Maximum of divergence of phiLumen after projection: "
                << max(mag(fvc::div(phiLumen_))) << endl
                << "Maximum of phiLumen: " << max(mag(phiLumen_)) << endl;
        }
    }
    
    // reconstruct cell-centered velocity
    UReconstructed_ = fvc::reconstruct(phiLumen_);

    if (debug)
    {
        if (time().outputTime())
        {
            UMappedLumen_.write();
            UReconstructed_.write();
            volScalarField divUReconstructed(fvc::div(phiLumen_));
            divUReconstructed.write();
        }
    }

    return phi_;
}


bool
Foam::graphVelocityEngine::writeData(Ostream& os) const
{
    os << "Lumen name: " << lumenMesh_.name() << nl
       << "Graph:      " << graph_ << endl;

    return os.good();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
