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

#include "plasmaRBCInletOutletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

#include "interpolation.H"

#include "RBCCollection.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

plasmaRBCInletOutletFvPatchScalarField::plasmaRBCInletOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchField<scalar>(p, iF),
    plasmaInletValue_(),
    RBCInletValue_(),
    plasmaFactor_(1.0),
    RBCFactor_(1.0),
    curTimeIndex_(-1)
{}


plasmaRBCInletOutletFvPatchScalarField::plasmaRBCInletOutletFvPatchScalarField
(
    const plasmaRBCInletOutletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchField<scalar>(ptf, p, iF, mapper),
    plasmaInletValue_(ptf.plasmaInletValue_),
    RBCInletValue_(ptf.RBCInletValue_),
    plasmaFactor_(ptf.plasmaFactor_),
    RBCFactor_(ptf.RBCFactor_),
    curTimeIndex_(-1)
{}


plasmaRBCInletOutletFvPatchScalarField::plasmaRBCInletOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchField<scalar>(p, iF, dict),
    plasmaInletValue_(readScalar(dict.lookup("plasmaInletValue"))),
    RBCInletValue_(readScalar(dict.lookup("RBCInletValue"))),
    plasmaFactor_(dict.lookupOrDefault<scalar>("plasmaFactor", 1.0)),
    RBCFactor_(dict.lookupOrDefault<scalar>("RBCFactor", 1.0)),
    curTimeIndex_(-1)
{}


plasmaRBCInletOutletFvPatchScalarField::plasmaRBCInletOutletFvPatchScalarField
(
    const plasmaRBCInletOutletFvPatchScalarField& ptf
)
:
    inletOutletFvPatchField<scalar>(ptf),
    plasmaInletValue_(ptf.plasmaInletValue_),
    RBCInletValue_(ptf.RBCInletValue_),
    plasmaFactor_(ptf.plasmaFactor_),
    RBCFactor_(ptf.RBCFactor_),
    curTimeIndex_(-1)
{}


plasmaRBCInletOutletFvPatchScalarField::plasmaRBCInletOutletFvPatchScalarField
(
    const plasmaRBCInletOutletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchField<scalar>(ptf, iF),
    plasmaInletValue_(ptf.plasmaInletValue_),
    RBCInletValue_(ptf.RBCInletValue_),
    plasmaFactor_(ptf.plasmaFactor_),
    RBCFactor_(ptf.RBCFactor_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void plasmaRBCInletOutletFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        const fvMesh& mesh = this->patch().boundaryMesh().mesh();
        const volScalarField& in_RBC_euler = mesh.lookupObject<volScalarField>("in_RBC_euler");
        const labelList& faceCells = this->patch().faceCells();

        // loop over all patch faces
        forAll(this->patch().Cf(), faceI)
        {
            label faceCellI = faceCells[faceI];

            if (in_RBC_euler[faceCellI] > 0.0)
            {
                this->refValue()[faceI] = RBCFactor_*RBCInletValue_;
            }
            else
            {
                this->refValue()[faceI] = plasmaFactor_*plasmaInletValue_;
            }
        }

        curTimeIndex_ = this->db().time().timeIndex();
    }

    inletOutletFvPatchField<scalar>::updateCoeffs();
}

void plasmaRBCInletOutletFvPatchScalarField::write(Ostream& os) const
{
    inletOutletFvPatchField<scalar>::write(os);
    os.writeKeyword("plasmaInletValue") << plasmaInletValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("RBCInletValue")    << RBCInletValue_    << token::END_STATEMENT << nl;
    os.writeKeyword("plasmaFactor")     << plasmaFactor_     << token::END_STATEMENT << nl;
    os.writeKeyword("RBCFactor")        << RBCFactor_        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    plasmaRBCInletOutletFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
