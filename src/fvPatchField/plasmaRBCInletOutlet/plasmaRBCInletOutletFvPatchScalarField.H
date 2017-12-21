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
    Foam::plasmaRBCInletOutletFvPatchScalarField

Description
    Modified inlet-outlet boundary condition that uses different inlet values
    when an RBC overlaps with the boundary.

SourceFiles
    plasmaRBCInletOutletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaRBCInletOutletFvPatchScalarField_H
#define plasmaRBCInletOutletFvPatchScalarField_H


#include "inletOutletFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class plasmaRBCInletOutletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class plasmaRBCInletOutletFvPatchScalarField
:
    public inletOutletFvPatchField<scalar>
{
    // Private data

        //- Field value at inlet when no RBC mesh overlaps it
        scalar plasmaInletValue_;

        //- Field value at inlet when a RBC mesh overlaps it
        scalar RBCInletValue_;

        //- Factor for values in plasma
        scalar plasmaFactor_;
        
        //- Factor for values in RBC
        scalar RBCFactor_;

        //- Current time index
        label curTimeIndex_;

public:

    //- Runtime type information
    TypeName("plasmaRBCInletOutlet");


    // Constructors

        //- Construct from patch and internal field
        plasmaRBCInletOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        plasmaRBCInletOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given plasmaRBCInletOutletFvPatchScalarField
        //  onto a new patch
        plasmaRBCInletOutletFvPatchScalarField
        (
            const plasmaRBCInletOutletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        plasmaRBCInletOutletFvPatchScalarField
        (
            const plasmaRBCInletOutletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new plasmaRBCInletOutletFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        plasmaRBCInletOutletFvPatchScalarField
        (
            const plasmaRBCInletOutletFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new plasmaRBCInletOutletFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return plasma PO2 at the inlet
            scalar plasmaInletValue() const
            {
                return plasmaInletValue_;
            }

            scalar& plasmaInletValue()
            {
                return plasmaInletValue_;
            }

            //- Return RBC PO2 at the inlet
            scalar RBCInletValue() const
            {
                return RBCInletValue_;
            }

            scalar& RBCInletValue()
            {
                return RBCInletValue_;
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

#endif

// ************************************************************************* //
