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

#include "analyticalSolutionPO2Cone.H"

#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

dimensionedScalar 
Foam::analyticalSolutionPO2Cone::clampedXCoordinate(const point& p) const
{
    scalar x = max(min(domainLength_.value(), p.x()), 0.0);
    return dimensionedScalar("dx", dimLength, x);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::analyticalSolutionPO2Cone::analyticalSolutionPO2Cone
(
    const dictionary& transportProperties,
    const dictionary& geometricProperties
)
:
    linearDensity_(transportProperties.lookup("linearDensity")),
    RBCVelocity_(transportProperties.lookup("RBCVelocity")),
    HbInlet_(transportProperties.lookup("HbInlet")),
    kappaO2RBC_(transportProperties.lookup("kappaO2RBC")),
    kappaO2Plasma_(transportProperties.lookup("kappaO2Plasma")),
    kappaO2Wall_(transportProperties.lookup("kappaO2Wall")),
    kappaO2Tissue_(transportProperties.lookup("kappaO2Tissue")),
    alphaRBC_(transportProperties.lookup("alphaRBC")),
    alphaPlasma_(transportProperties.lookup("alphaPlasma")),
    alphaWall_(transportProperties.lookup("alphaWall")),
    alphaTissue_(transportProperties.lookup("alphaTissue")),
    VMolO2_(transportProperties.lookup("VMolO2")),
    NHb_(transportProperties.lookup("NHb")),
    O2ConsumptionRate_(transportProperties.lookup("O2ConsumptionRateConstant")),
    radiusRBC_(geometricProperties.lookup("radiusRBC")),
    radiusPlasma_(geometricProperties.lookup("radiusPlasma")),
    radiusWall_(geometricProperties.lookup("radiusWall")),
    radiusTissueLeft_(geometricProperties.lookup("radiusTissueLeft")),
    radiusTissueRight_(geometricProperties.lookup("radiusTissueRight")),
    domainLength_(geometricProperties.lookup("domainLength")),
    RBCLength_(geometricProperties.lookup("RBCLength")),
    RBCVolume_(geometricProperties.lookup("RBCVolume")),
    dissociationCurvePtr_
    (
        dissociationCurve::New(transportProperties.subDict("dissociationCurve"))
    ),
    minimalSaturation_
    (
        transportProperties.lookupOrDefault
        (
            "minimalSaturation",
            dimensionedScalar("minimalSaturation", dimless, 0.01)
        )
    ),
    minimalPO2_
    (
        transportProperties.lookupOrDefault
        (
            "minimalPO2",
            dimensionedScalar("minimalPO2", dimless, 1.)
        )
    ),
    intravascResistanceCoeffLDHalf_
    (
        transportProperties.lookupOrDefault
        (
            "IVRLDHalf", 
            dimensionedScalar("IVRLDHalf", dimLength*dimTime, 0.0)
        )
    ),
    useAnalyticalIntravascResistanceCoeff_
    (
        readBool(transportProperties.lookup("useAnalyticalIVR"))
    )
{
    if (!useAnalyticalIntravascResistanceCoeff_ &&
        !transportProperties.found("IVRLDHalf"))
    {
        FatalErrorIn
        (
            "analyticalSolutionPO2Cone::analyticalSolutionPO2Cone("
            "const dictionary&, const dictionary&)"
        )   << "  A prescribed IVR coefficient should be given but none is given." << nl
            << "  A coefficient should be given under the key ""IVRLDHalf""." << nl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::analyticalSolutionPO2Cone::~analyticalSolutionPO2Cone()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

dimensionedScalar 
Foam::analyticalSolutionPO2Cone::intravascularResistanceCoefficient() const
{
    if (useAnalyticalIntravascResistanceCoeff_)
    {
        return 1./(2.*constant::mathematical::pi*linearDensity_)
               * (   1./(4.*kappaO2RBC_*alphaRBC_)
                   + log(radiusPlasma_/radiusRBC_)   /(kappaO2Plasma_*alphaPlasma_)
                   + log(radiusWall_  /radiusPlasma_)/(kappaO2Wall_  *alphaWall_  )
                 );
    }
    else
    {
        return 0.5/linearDensity_ * intravascResistanceCoeffLDHalf_;
    }
}


dimensionedScalar
Foam::analyticalSolutionPO2Cone::RBCFlux() const
{
    return linearDensity_*RBCVelocity_/RBCLength_;
}


dimensionedScalar
Foam::analyticalSolutionPO2Cone::radiusAtX(const dimensionedScalar& dx) const
{
    return (1. - dx/domainLength_)*radiusTissueLeft_ + dx/domainLength_*radiusTissueRight_;
}


dimensionedScalar
Foam::analyticalSolutionPO2Cone::consumptionToX(const dimensionedScalar& dx) const
{
    dimensionedScalar R = radiusAtX(dx);
    return dx*constant::mathematical::pi*O2ConsumptionRate_
           * (1./3.*(sqr(radiusTissueLeft_) + radiusTissueLeft_*R + sqr(R)) - sqr(radiusWall_));
}


dimensionedScalar
Foam::analyticalSolutionPO2Cone::consumptionPerLengthAtX(const dimensionedScalar& dx) const
{
    return constant::mathematical::pi*(sqr(radiusAtX(dx)) - sqr(radiusWall_))*O2ConsumptionRate_;
}


dimensionedScalar
Foam::analyticalSolutionPO2Cone::saturationAtX(const dimensionedScalar& dx) const
{
    dimensionedScalar convO2Capacity = RBCFlux()*VMolO2_*NHb_*RBCVolume_;
    dimensionedScalar S = HbInlet_ - consumptionToX(dx)/convO2Capacity;
    S.value() = max(S.value(), minimalSaturation_.value());
    return S;
}


dimensionedScalar
Foam::analyticalSolutionPO2Cone::convectivePO2Drop(const dimensionedScalar& dx) const
{
    return dissociationCurvePtr_->equilibriumPO2(HbInlet_) 
         - dissociationCurvePtr_->equilibriumPO2(saturationAtX(dx));
}


dimensionedScalar
Foam::analyticalSolutionPO2Cone::intravascResistancePO2Drop(const dimensionedScalar& dx) const
{
    return intravascularResistanceCoefficient()*consumptionPerLengthAtX(dx);
}


dimensionedScalar
Foam::analyticalSolutionPO2Cone::extravascularPO2Drop(const point& p) const
{
    dimensionedScalar dx = clampedXCoordinate(p);
    dimensionedScalar dy("dy", dimLength, p.y());
    dimensionedScalar tissueRadius = radiusAtX(dx);
    return -O2ConsumptionRate_/(4.*kappaO2Tissue_*alphaTissue_)
           * (sqr(dy) - sqr(radiusWall_) - 2*sqr(tissueRadius)*log(dy/radiusWall_));
}


dimensionedScalar
Foam::analyticalSolutionPO2Cone::PO2Plasma(const point& p) const
{
    dimensionedScalar dx = clampedXCoordinate(p);
    return dissociationCurvePtr_->equilibriumPO2(saturationAtX(dx));
}


dimensionedScalar
Foam::analyticalSolutionPO2Cone::PO2Wall(const point& p) const
{
    dimensionedScalar dx = clampedXCoordinate(p);
    dimensionedScalar PO2Wall = dissociationCurvePtr_->equilibriumPO2(HbInlet_) 
                              - convectivePO2Drop(dx)
                              - intravascResistancePO2Drop(dx);
    PO2Wall.value() = max(PO2Wall.value(), minimalPO2_.value());
    return PO2Wall;
}


dimensionedScalar
Foam::analyticalSolutionPO2Cone::PO2Tissue(const point& p) const
{
    dimensionedScalar PO2Tissue = PO2Wall(p) - extravascularPO2Drop(p);
    PO2Tissue.value() = max(PO2Tissue.value(), minimalPO2_.value());
    return PO2Tissue;
}


void Foam::analyticalSolutionPO2Cone::computePO2
(
    volScalarField& PO2
) const
{
    const fvMesh& mesh = PO2.mesh();
    forAll(mesh.C(), cI)
    {
        const point p = mesh.C()[cI];
        if (p.y() >= radiusWall_.value())
        {
            PO2[cI] = PO2Tissue(p).value();
        }
        else if (p.y() >= radiusPlasma_.value())
        {
            PO2[cI] = PO2Wall(p).value();
        }
        else
        {
            PO2[cI] = PO2Plasma(p).value();
        }
    }
}


void Foam::analyticalSolutionPO2Cone::computeHbAndRBCPO2
(
    RBCCollection& RBCs
) const
{
    const fvMesh& RBCMesh = RBCs.RBCMesh();
    volScalarField& Hb = RBCs.Hb();
    volScalarField& PO2RBC = RBCs.PO2();
    forAll(RBCMesh.C(), cI)
    {
        point p = RBCMesh.C()[cI];
        dimensionedScalar dx = clampedXCoordinate(p);
        Hb[cI]  = saturationAtX(dx).value();
        PO2RBC[cI] = dissociationCurvePtr_->equilibriumPO2(Hb[cI]);
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::analyticalSolutionPO2Cone::operator=(const analyticalSolutionPO2Cone& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::analyticalSolutionPO2Cone::operator=(const Foam::analyticalSolutionPO2Cone&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// ************************************************************************* //
