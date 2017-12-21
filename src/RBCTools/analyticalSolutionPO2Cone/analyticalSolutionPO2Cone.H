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
    Foam::analyticalSolutionPO2Cone

Description

SourceFiles
    analyticalSolutionPO2Cone.C

\*---------------------------------------------------------------------------*/

#ifndef analyticalSolutionPO2Cone_H
#define analyticalSolutionPO2Cone_H

#include "fvCFD.H"

#include "RBCCollection.H"
#include "dissociationCurve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class analyticalSolutionPO2Cone Declaration
\*---------------------------------------------------------------------------*/

class analyticalSolutionPO2Cone
{
    // Private data

        dimensionedScalar linearDensity_;
        dimensionedScalar RBCVelocity_;

        dimensionedScalar HbInlet_;
        dimensionedScalar kappaO2RBC_;
        dimensionedScalar kappaO2Plasma_;
        dimensionedScalar kappaO2Wall_;
        dimensionedScalar kappaO2Tissue_;
        dimensionedScalar alphaRBC_;
        dimensionedScalar alphaPlasma_;
        dimensionedScalar alphaWall_;
        dimensionedScalar alphaTissue_;

        dimensionedScalar VMolO2_;
        dimensionedScalar NHb_;

        dimensionedScalar O2ConsumptionRate_;

        dimensionedScalar radiusRBC_;
        dimensionedScalar radiusPlasma_;
        dimensionedScalar radiusWall_;
        dimensionedScalar radiusTissueLeft_;
        dimensionedScalar radiusTissueRight_;
        dimensionedScalar domainLength_;
        dimensionedScalar RBCLength_;
        dimensionedScalar RBCVolume_;

        autoPtr<dissociationCurve> dissociationCurvePtr_;

        dimensionedScalar minimalSaturation_;
        dimensionedScalar minimalPO2_;
        dimensionedScalar intravascResistanceCoeffLDHalf_;
        bool useAnalyticalIntravascResistanceCoeff_;

    // Private Member Functions

        //- Return the clamped x-coordinate so that 0 <= x <= domainLength_.
        dimensionedScalar clampedXCoordinate(const point& p) const;

        //- Disallow default bitwise copy construct
        analyticalSolutionPO2Cone(const analyticalSolutionPO2Cone&);

        //- Disallow default bitwise assignment
        void operator=(const analyticalSolutionPO2Cone&);


public:

    // Constructors

        //- Construct from dictionaries
        analyticalSolutionPO2Cone
        (
            const dictionary& oxygenProperties,
            const dictionary& geometricProperties
        );


    //- Destructor
    ~analyticalSolutionPO2Cone();


    // Member Functions

        //- Compute the analytical intravascular resistance coefficient
        dimensionedScalar intravascularResistanceCoefficient() const;

        dimensionedScalar RBCFlux() const;
        dimensionedScalar radiusAtX(const dimensionedScalar& dx) const;
        dimensionedScalar consumptionToX(const dimensionedScalar& dx) const;
        dimensionedScalar consumptionPerLengthAtX(const dimensionedScalar& dx) const;
        dimensionedScalar saturationAtX(const dimensionedScalar& dx) const;
        dimensionedScalar convectivePO2Drop(const dimensionedScalar& dx) const;
        dimensionedScalar intravascResistancePO2Drop(const dimensionedScalar& dx) const;
        dimensionedScalar extravascularPO2Drop(const point& p) const;
        dimensionedScalar PO2Plasma(const point& p) const;
        dimensionedScalar PO2Wall(const point& p) const;
        dimensionedScalar PO2Tissue(const point& p) const;

        //- Compute the analytical solution for PO2
        void computePO2(volScalarField& PO2) const;

        //- Compute the analytical solution for hemoglobin saturation
        void computeHbAndRBCPO2(RBCCollection&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
