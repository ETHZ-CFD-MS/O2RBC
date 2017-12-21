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

#include "RBCPathMover.H"

#include "RBC.H"
#include "RBCPath.H"
#include "RBCCollection.H"
#include "RBCPathCollection.H"

#include "circularTubeGraph.H"
#include "axisymmetricBodyGeometricState.H"
#include "diameterFunction.H"
#include "dissociationCurve.H"
#include "graphInletValue.H"

#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RBCPathMover, 0);

    addToRunTimeSelectionTable
    (
        RBCMover,
        RBCPathMover,
        meshRBCCollection
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

autoPtr<deformableBodyGeometricState>
RBCPathMover::computeRBCGeometricState(const RBCPath& path) const
{
    // interpolate the graph coordinate
    graphCoordinateDirected gc = graphCoordinate(path.index());
    label eI          = gc.edgeIndex();
    scalar s          = gc.sCoord();
    bool goingForward = gc.goingForward();

    // set the RBC position
    point center = graph_.edge(eI).pointPosition(s);

    // set RBC direction
    vector axis = graph_.edge(eI).tangentVector(s);
    if (!goingForward)
    {
        axis = -axis;
    }

    scalar diameter = RBCDiameterFunctionPtr_()(graph_.tube(eI).innerDiameter(s));
    scalar length   = RBCVolume_/(constant::mathematical::pi*sqr(diameter/2.0));

    return autoPtr<deformableBodyGeometricState>
        (
            new axisymmetricBodyGeometricState(center, axis, diameter, length)
        );
}


scalar
RBCPathMover::computeInletValue(const RBCPath& path) const
{
    const label eI = path.edges()[0];
    const scalar inletValue = 
        inletValuePtr_->inletValue(eI, mesh_.time().timeOutputValue());
    const word inletFieldName = inletValuePtr_->fieldName();
    return computePO2FromFieldNameAndValue(inletValue, inletFieldName);
}


scalar 
RBCPathMover::currentInletValue(const RBCPath& path) const
{
    graphCoordinateDirected gc = graphInterpolator_.interpolate
    (
        path.pathTimes(),
        path.edges(),
        path.sCoords(),
        mesh_.time().timeOutputValue()
    );
    return inletValuePtr_->inletValue(gc.edgeIndex(), mesh_.time().timeOutputValue());
}


void
RBCPathMover::normalizePathCoordinates()
{
    forAll(RBCPaths_, i)
    {
        RBCPath& cPath = RBCPaths_[i];

        forAll(cPath.sCoords(), cI)
        {
            scalar edgeLength = graph_.edge(cPath.edges()[cI]).length();
            cPath.scaleSCoord(cI, 1.0/edgeLength);
        }
        if (debug)
        {
            Info<< "Curvilinear coordinates of path " << cPath.index() << ": "
                << cPath.sCoords() << endl;
        }
    }
}


void
RBCPathMover::convertPathEdgeIndices()
{
    if (debug)
    {
        Info<< "RBCPathMover::convertPathEdgeIndices(): "
            << "Converting external to internal path edge indices." << endl;
    }

    forAll(RBCPaths_, i)
    {
        RBCPath& cPath = RBCPaths_[i];

        forAll(cPath.edges(), i)
        {
            label externalEdgeI = cPath.edges()[i];
            label internalEdgeI = graph_.externalToInternalEdgeIndex(externalEdgeI);

            cPath.setEdge(i, internalEdgeI);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RBCPathMover::RBCPathMover
(
    const fvMesh& mesh,
    RBCCollection& RBCCollection
)
:
    RBCMover(mesh, RBCCollection),
    RBCPaths_
    (
        IOobject
        (
            "RBCPaths",
            mesh.time().caseConstant(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    graph_(mesh_.time().lookupObject<circularTubeGraph>("graphDict")),
    graphInterpolator_(graph_),
    sampleMeshNames_(subDict("sampleMeshNames")),
    RBCDiameterFunctionPtr_
    (
        diameterFunction::New(subDict("diameterFunction"))
    ),
    inletValuePtr_
    (
        graphInletValue:: New(subDict("inletValues"), graph_)
    ),
    initialFieldName_(lookup("initialField")),
    initialValue_(readScalar(lookup("initialValue"))),
    RBCVolume_(readScalar(lookup("RBCVolume")))
{

    if (lookupOrDefault<bool>("useExternalEdgeIndices", true))
    {
        Info<< "Converting external to internal edge indices..." << endl;
        convertPathEdgeIndices();
    }

    if (readBool(lookup("normalizePathCoordinates")))
    {
        Info<< "Normalizing path coordinates..." << endl;
        normalizePathCoordinates();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

RBCPathMover::~RBCPathMover()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
RBCPathMover::setInitialPositions()
{
    labelList activePathIdx = RBCPaths_.activeIndices();

    if (debug)
    {
        Info<< "Starting to set initial positions" << endl;
        Info<< "Number of RBC paths: " << RBCPaths_.size() << endl;
        Info<< "Number of active RBC paths: " << activePathIdx.size() << endl;
    }

    // attach paths to RBCs
    forAll(activePathIdx, i)
    {
        label pathI = activePathIdx[i];
        autoPtr<deformableBodyGeometricState> geometricStatePtr
                            = computeRBCGeometricState(RBCPaths_(pathI));
        label edgeI = graphCoordinate(pathI).edgeIndex();
        word sampleMeshName = sampleMeshNames_.value(pathI, edgeI);
        RBCCollection_.insert
        (
            RBC
            (
                pathI, 
                geometricStatePtr(), 
                sampleMeshName,
                computeInletValue(RBCPaths_(pathI))
            )
        );
    }
    RBCCollection_.updateMeshes();
}


void
RBCPathMover::moveAll()
{
    // remove RBCs corresponding to new inactive paths
    labelList newInactivePaths = RBCPaths_.newInactiveIndices();

    forAll(newInactivePaths, i)
    {
        label pathI = newInactivePaths[i];
        RBCCollection_.remove(pathI);
    }

    // insert RBCs for new active paths
    labelList newActivePaths = RBCPaths_.newActiveIndices();
    
    forAll(newActivePaths, i)
    {
        label pathI = newActivePaths[i];
        autoPtr<deformableBodyGeometricState> geometricStatePtr
                            = computeRBCGeometricState(RBCPaths_(pathI));
        label edgeI = graphCoordinate(pathI).edgeIndex();
        word sampleMeshName = sampleMeshNames_.value(pathI, edgeI);
        RBCCollection_.insert
        (
            RBC
            (
                pathI, 
                geometricStatePtr(), 
                sampleMeshName,
                computeInletValue(RBCPaths_(pathI))
            )
        );
    }

    // move all RBCs
    labelList activePaths = RBCPaths_.activeIndices();
    forAll(activePaths, i)
    {
        label pathI = activePaths[i];
        autoPtr<deformableBodyGeometricState> geometricStatePtr
                            = computeRBCGeometricState(RBCPaths_(pathI));
        RBCCollection_[pathI].setInletPO2
        (
            RBCCollection_.PO2FromFieldNameAndValue
            (
                inletValuePtr_->fieldName(),
                currentInletValue(RBCPaths_(pathI))
            )
        );
        RBCCollection_.prepareMotion(pathI, geometricStatePtr());
    }
    RBCCollection_.applyMotion();
}


graphCoordinateDirected
Foam::RBCPathMover::graphCoordinate(const label pathI) const
{
    const RBCPath& path = RBCPaths_(pathI);
    return graphInterpolator_.interpolate
    (
        path.pathTimes(),
        path.edges(),
        path.sCoords(),
        mesh_.time().timeOutputValue()
    );
}


Foam::Ostream&
RBCPathMover::write
(
    Ostream& os
) const
{
    os.writeKeyword("sampleMeshNames")  << nl << sampleMeshNames_ << nl;
    os.writeKeyword("inletValues")      << nl << inletValuePtr_() << nl;
    os.writeKeyword("diameterFunction") << nl << RBCDiameterFunctionPtr_() << nl;
    os.writeKeyword("initialField") << initialFieldName_<< token::END_STATEMENT << nl;
    os.writeKeyword("initialValue") << initialValue_    << token::END_STATEMENT << nl;
    os.writeKeyword("RBCVolume")    << RBCVolume_       << token::END_STATEMENT << nl;

    return os;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream& 
Foam::operator<<(Ostream& os, const RBCPathMover& obj)
{
    obj.write(os);

    return os;
}

// ************************************************************************* //
