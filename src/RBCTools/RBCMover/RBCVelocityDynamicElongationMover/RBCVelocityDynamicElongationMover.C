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

#include "RBCVelocityDynamicElongationMover.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RBCVelocityDynamicElongationMover, 0);
    const word RBCVelocityDynamicElongationMover::plasmaZoneName = "plasma";

    addToRunTimeSelectionTable
    (
        RBCMover,
        RBCVelocityDynamicElongationMover,
        meshRBCCollection
    );

    addToRunTimeSelectionTable
    (
        RBCVelocityMover,
        RBCVelocityDynamicElongationMover,
        meshRBCCollection
    );
}

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

Foam::scalar RBCVelocityDynamicElongationMover::PO2Average() const
{
    const volScalarField& inRBCEuler = mesh_.lookupObject<volScalarField>("inRBCEuler");
    const volScalarField& PO2 = mesh_.lookupObject<volScalarField>("PO2");
    const cellZoneMesh& czMesh = mesh_.cellZones();
    const cellZone& plasmaZone = czMesh[RBCVelocityDynamicElongationMover::plasmaZoneName];
    const volScalarField::Internal& V = mesh_.V();
    scalar plasmaVolume = 0.0;
    scalar PO2SumPlasma = 0.0;
    forAll(plasmaZone, cI)
    {
        plasmaVolume += (1.0 - inRBCEuler[cI])*V[cI];
        PO2SumPlasma += PO2[cI]*(1.0 - inRBCEuler[cI])*V[cI];
    }
    return PO2SumPlasma/plasmaVolume;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


RBCVelocityDynamicElongationMover::RBCVelocityDynamicElongationMover
(
    const fvMesh& mesh,
    RBCCollection& RBCCollection
)
:
    RBCVelocityMover(mesh, RBCCollection),
    tubeHematocrit_(),
    RBCVelocityReference_(RBCVelocity_),
    RBCRadiusTable_(subDict("RBCRadiusTable")),
    RBCVelocityFactorTable_(subDict("RBCVelocityFactor")),
    PO2AveragingZone_(lookup("PO2AveragingZone"))
{
    tubeHematocrit_ = linearDensity_*sqr(0.5*RBCDiameter()/plasmaRadius_);
}


RBCVelocityDynamicElongationMover::~RBCVelocityDynamicElongationMover()
{}


void RBCVelocityDynamicElongationMover::moveAll()
{
    RBCVelocity_ = RBCVelocityFactorTable_(PO2Average())*RBCVelocityReference_;
    linearDensity_ = tubeHematocrit_*sqr(plasmaRadius_/RBCRadius());
    if (debug)
    {
        Info<< "Foam::RBCVelocityDynamicElongationMover::moveAll():" << endl;
        Info<< "    Average PO2 in the plasma: " << PO2Average() << endl;
        Info<< "    RBC velocity:   " << RBCVelocity_ << endl;
        Info<< "    Linear density: " << linearDensity_ << endl;
        Info<< "    RBC diameter:   " << RBCDiameter() << endl;
    }
    RBCVelocityMover::moveAll();
}


scalar RBCVelocityDynamicElongationMover::RBCRadius() const
{
    return RBCRadiusTable_(PO2Average());
}

