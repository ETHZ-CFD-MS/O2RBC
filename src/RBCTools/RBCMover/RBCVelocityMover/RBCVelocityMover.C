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

#include "RBCVelocityMover.H"

#include "RBCCollection.H"

#include "axisymmetricBodyGeometricState.H"
#include "randomVariable.H"

#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RBCVelocityMover, 0);
    defineRunTimeSelectionTable(RBCVelocityMover, meshRBCCollection);

    addToRunTimeSelectionTable
    (
        RBCMover,
        RBCVelocityMover,
        meshRBCCollection
    );
    addToRunTimeSelectionTable
    (
        RBCVelocityMover,
        RBCVelocityMover,
        meshRBCCollection
    );
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::RBCVelocityMover> Foam::RBCVelocityMover::New
(
    const fvMesh& mesh,
    RBCCollection& RBCCollection
)
{
    const word functionType
    (
        IOdictionary
        (
            IOobject
            (
                "RBCMoverDict",
                mesh.time().system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("type")
    );

    meshRBCCollectionConstructorTable::iterator cstrIter =
        meshRBCCollectionConstructorTablePtr_->find(functionType);

    if (cstrIter == meshRBCCollectionConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "RBCVelocityMover::New(const fvMesh&, const RBCCollection&)"
        )   << "Unknown RBCVelocityMover type " << functionType
            << nl << nl
            << "Valid RBCVelocityMover types : " << endl
            << meshRBCCollectionConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<RBCVelocityMover>(cstrIter()(mesh, RBCCollection));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


scalar
RBCVelocityMover::computeInletValue() const
{
    const scalar inletValue = inletValuePtr_->generate();
    return computePO2FromFieldNameAndValue(inletValue, inletFieldName_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

RBCVelocityMover::RBCVelocityMover
(
    const fvMesh& mesh,
    RBCCollection& RBCCollection
)
:
    RBCMover(mesh, RBCCollection),
    inletValuePtr_
    (
        new randomVariable(subDict("inletValue"))
    ),
    inletFieldName_(lookup("inletFieldName")),
    linearDensity_(readScalar(lookup("linearDensity"))),
    RBCVelocity_(readScalar(lookup("RBCVelocity"))),
    initialFieldName_(lookup("initialField")),
    initialValue_(readScalar(lookup("initialValue"))),
    RBCRadius_(readScalar(lookup("RBCRadius"))),
    RBCVolume_(readScalar(lookup("RBCVolume"))),
    plasmaRadius_(readScalar(lookup("plasmaRadius"))),
    sampleMeshName_(lookup("sampleMeshName"))
{}


RBCVelocityMover::~RBCVelocityMover()
{}


void RBCVelocityMover::setInitialPositions()
{
    const vector axis(1.0, 0.0, 0.0);
    const scalar xMin = mesh_.bounds().min().x();
    const scalar xMax = mesh_.bounds().max().x();
    scalar RBCLength = RBCVolume_/(constant::mathematical::pi*sqr(RBCRadius()));

    // The RBCs are numbered from right to left, since the RBC index should
    // correspond to the order of the entry times into the domain.
    // Therefore, the number of RBCs in the domain is first determined.
    scalar xRBC = xMin;
    label nRBC = 0;
    while (xRBC < xMax + 0.5*RBCLength)
    {
        xRBC += RBCLength/linearDensity_;
        nRBC++;
    }

    xRBC = xMin;
    for (int RBCI = nRBC - 1; RBCI >= 0; RBCI--)
    {
        point center(xRBC, 0.0, 0.0);
        autoPtr<deformableBodyGeometricState> geometricStatePtr
        (
            new axisymmetricBodyGeometricState(center, axis, RBCDiameter(), RBCLength)
        );
        RBCCollection_.insert
        (
            RBC
            (
                RBCI, 
                geometricStatePtr(), 
                sampleMeshName_,
                computePO2FromFieldNameAndValue(initialValue_, initialFieldName_)
            )
        );
        xRBC += RBCLength/linearDensity_;
    }
    RBCCollection_.updateMeshes();
}


void RBCVelocityMover::moveAll()
{
    const vector axis(1.0, 0.0, 0.0);
    const scalar dt = mesh_.time().deltaT().value();
    scalar RBCLength = RBCVolume_/(constant::mathematical::pi*sqr(RBCRadius()));
    scalar minimalXCenter = VGREAT;
    label maximalRBCI = 0;
    labelList removeList;

    forAll(RBCCollection_.RBCIndices(), i)
    {
        // advance existing RBCs and remove them they left the domain
        label RBCI = RBCCollection_.RBCIndices()[i];
        RBC& cRBC = RBCCollection_[RBCI];
        const deformableBodyGeometricState& oldGeometricState = cRBC.geometricState();
        scalar newXCenter = oldGeometricState.center().x() + dt*RBCVelocity_;

        if (newXCenter < mesh_.bounds().max().x() + 0.5*RBCLength
         && newXCenter >= mesh_.bounds().min().x() - 0.5*RBCLength)
        {
            autoPtr<deformableBodyGeometricState> geometricStatePtr
            (
                new axisymmetricBodyGeometricState
                (
                    point(newXCenter, 0.0, 0.0), 
                    axis, 
                    RBCDiameter(), 
                    RBCLength
                )
            );
            RBCCollection_.prepareMotion(cRBC.index(), geometricStatePtr());
            minimalXCenter = min(minimalXCenter, newXCenter);
            maximalRBCI = max(maximalRBCI, cRBC.index());
        }
        else
        {
            removeList.append(cRBC.index());
        }
    }
    // remove RBCs 
    forAll(removeList, i)
    {
        RBCCollection_.remove(removeList[i]);
    }

    // add new RBC if needed
    scalar RBCSpacing = RBCLength/linearDensity_;
    if (minimalXCenter + 0.5*RBCLength >= mesh_.bounds().min().x() + RBCSpacing)
    {
        autoPtr<deformableBodyGeometricState> geometricStatePtr
        (
            new axisymmetricBodyGeometricState
            (
                point(minimalXCenter - RBCSpacing, 0.0, 0.0), 
                axis, 
                RBCDiameter(), 
                RBCLength
            )
        );
        RBCCollection_.insert
        (
            RBC
            (
                maximalRBCI + 1, 
                geometricStatePtr(), 
                sampleMeshName_,
                computeInletValue()
            )
        );
    }
    RBCCollection_.applyMotion();
}


scalar RBCVelocityMover::RBCDiameter() const
{
    return 2.0*RBCRadius();
}


scalar RBCVelocityMover::RBCRadius() const
{
    return RBCRadius_;
}


bool RBCVelocityMover::writeData
(
    Ostream& os
) const
{
    return os.good();
}

