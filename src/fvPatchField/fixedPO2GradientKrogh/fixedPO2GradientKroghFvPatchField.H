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
    Foam::fixedPO2GradientKroghFvPatchField

Description
    Fixed gradient boundary condition for PO2 based on the Krogh cylinder
    solution to the oxygen transport equation. This class uses a moving frame
    for the RBC, as the inlet hemoglobin saturation is a function of time 
    and RBC velocity.

    A sample entry in a field's boundaryField is given by

    patchName
    {
        type                fixedPO2GradientKrogh;
        P50                 47.9;    // equilibrium PO2 at half saturation
        hillExponent        2.5;     // Hill exponent
        S_inlet             0.8;     // hemoglobin saturation at the inlet
        z_inlet             0;       // axial coordinate of the inlet
        RBC_velocity        1e-3;    // RBC velocity
        R_c                 1.5e-6;  // RBC radius
        R_w                 2.6e-6;  // capillary wall radius
        R_t                 20e-6;   // tissue radius
        O2ConsumptionRate   1000;    // metabolic rate of oxygen consumption
        hematocrit          0.3;     // tube hematocrit
        NHb                 20.3;    // heme group concentration
        VMolO2              25440;   // molar volume of oxygen
    }

SourceFiles
    fixedPO2GradientKroghFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedPO2GradientKroghFvPatchField_H
#define fixedPO2GradientKroghFvPatchField_H


#include "fixedGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class fixedPO2GradientKroghFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedPO2GradientKroghFvPatchField
:
    public fixedGradientFvPatchField<Type>
{
    // Private data

        //- PO2 at 50% saturation
        Type P50_;

        //- Hemoglobin saturation at inlet
        scalar S_inlet_;

        //- axial coordinate of the inlet
        scalar z_inlet_;

        //- Coefficient in axial saturation equation
        scalar beta_;

        //- Hill exponent
        dimensionedScalar hillExponent_;

        //- Moving frame velocity (= RBC velocity)
        dimensionedScalar RBC_velocity_;

        //- Current time index
        label curTimeIndex_;

public:

    //- Runtime type information
    TypeName("fixedPO2GradientKrogh");


    // Constructors

        //- Construct from patch and internal field
        fixedPO2GradientKroghFvPatchField<Type>
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedPO2GradientKroghFvPatchField<Type>
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given fixedPO2GradientKroghFvPatchField
        //  onto a new patch
        fixedPO2GradientKroghFvPatchField<Type>
        (
            const fixedPO2GradientKroghFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedPO2GradientKroghFvPatchField<Type>
        (
            const fixedPO2GradientKroghFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new fixedPO2GradientKroghFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedPO2GradientKroghFvPatchField<Type>
        (
            const fixedPO2GradientKroghFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new fixedPO2GradientKroghFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return inlet hemoglobin saturation
            scalar S_inlet() const
            {
                return S_inlet_;
            }

            scalar& S_inlet()
            {
                return S_inlet_;
            }

            //- Return inlet axial coordinate
            scalar z_inlet() const
            {
                return z_inlet_;
            }

            scalar& z_inlet()
            {
                return z_inlet_;
            }

            //- Return beta
            scalar beta() const
            {
                return beta_;
            }

            scalar& beta()
            {
                return beta_;
            }

            //- Return P50
            Type P50() const
            {
                return P50_;
            }

            Type& P50()
            {
                return P50_;
            }

            //- Return hillExponent
            dimensionedScalar hillExponent() const
            {
                return hillExponent_;
            }

            dimensionedScalar& hillExponent()
            {
                return hillExponent_;
            }

            //- Return v
            dimensionedScalar RBC_velocity() const
            {
                return RBC_velocity_;
            }

            dimensionedScalar& RBC_velocity()
            {
                return RBC_velocity_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            
        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fixedPO2GradientKroghFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
