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

#include "fixedPO2GradientKroghFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
fixedPO2GradientKroghFvPatchField<Type>::fixedPO2GradientKroghFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    P50_(pTraits<Type>::zero),
    S_inlet_(),
    z_inlet_(),
    beta_(),
    hillExponent_(),
    RBC_velocity_(),
    curTimeIndex_(-1)
{}


template<class Type>
fixedPO2GradientKroghFvPatchField<Type>::fixedPO2GradientKroghFvPatchField
(
    const fixedPO2GradientKroghFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<Type>(ptf, p, iF, mapper),
    P50_(ptf.P50_),
    S_inlet_(ptf.S_inlet_),
    z_inlet_(ptf.z_inlet_),
    beta_(ptf.beta_),
    hillExponent_(ptf.hillExponent_),
    RBC_velocity_(ptf.RBC_velocity_),
    curTimeIndex_(-1)
{}


template<class Type>
fixedPO2GradientKroghFvPatchField<Type>::fixedPO2GradientKroghFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<Type>(p, iF, dict),
    P50_(dict.lookupOrDefault<Type>("P50", pTraits<Type>::zero)),
    S_inlet_(dict.lookupOrDefault<scalar>("S_inlet", 0.8)),
    z_inlet_(dict.lookupOrDefault<scalar>("z_inlet",-2e-6)),
    beta_(),
    hillExponent_("hillExponent", dimensionSet(0,0,0,0,0,0,0), 2.2),
    RBC_velocity_("RBC_velocity", dimensionSet(0,1,-1,0,0,0,0), 0.0),
    curTimeIndex_(-1)
{
    hillExponent_ = dict.lookup("hillExponent");
    RBC_velocity_ = dict.lookup("RBC_velocity");

    if (dict.found("beta"))
    {
        beta_ = dict.lookupOrDefault<scalar>("beta", 500.0);
    }
    else
    {
        // compute beta from values in dictionary
        dimensionedScalar R_c("R_c", dimensionSet(0,1,0,0,0,0,0), dict.lookup("R_c"));
        dimensionedScalar R_w("R_t", dimensionSet(0,1,0,0,0,0,0), dict.lookup("R_w"));
        dimensionedScalar R_t("R_w", dimensionSet(0,1,0,0,0,0,0), dict.lookup("R_t"));
        dimensionedScalar O2ConsumptionRate
        (
            "O2ConsumptionRate", 
            dimensionSet(0,0,-1,0,0,0,0),
            dict.lookup("O2ConsumptionRate")
        );
        dimensionedScalar hematocrit("hematocrit", dimensionSet(0,0,0,0,0,0,0), dict.lookup("hematocrit"));
        dimensionedScalar NHb("NHb", dimensionSet(0,-3,0,0,0,0,0),dict.lookup("NHb"));
        dimensionedScalar VMolO2("VMolO2", dimensionSet(0,0,0,0,0,0,0), dict.lookup("VMolO2"));

        dimensionedScalar beta_dim = O2ConsumptionRate*(R_t*R_t - R_w*R_w)
                           / (hematocrit*NHb*VMolO2*RBC_velocity_*R_c*R_c);

        beta_ = beta_dim.value();
    }
}


template<class Type>
fixedPO2GradientKroghFvPatchField<Type>::fixedPO2GradientKroghFvPatchField
(
    const fixedPO2GradientKroghFvPatchField& ptf
)
:
    fixedGradientFvPatchField<Type>(ptf),
    P50_(ptf.P50_),
    S_inlet_(ptf.S_inlet_),
    z_inlet_(ptf.z_inlet_),
    beta_(ptf.beta_),
    hillExponent_(ptf.hillExponent_),
    RBC_velocity_(ptf.RBC_velocity_),
    curTimeIndex_(-1)
{}


template<class Type>
fixedPO2GradientKroghFvPatchField<Type>::fixedPO2GradientKroghFvPatchField
(
    const fixedPO2GradientKroghFvPatchField& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(ptf, iF),
    P50_(ptf.P50_),
    S_inlet_(ptf.S_inlet_),
    z_inlet_(ptf.z_inlet_),
    beta_(ptf.beta_),
    hillExponent_(ptf.hillExponent_),
    RBC_velocity_(ptf.RBC_velocity_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void fixedPO2GradientKroghFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<Type>& gradient = this->gradient();

        // compute axial coordinate of patch
        boundBox bb(this->patch().patch().localPoints(), true);
        const vector ctr = 0.5*(bb.max() + bb.min());
        const scalar   z = ctr & vector(1,0,0);

        const scalar t = this->db().time().timeOutputValue();

        // hemoglobin saturation at (z,t) (moving frame)
        const scalar S = S_inlet_ - beta_*(z - z_inlet_ + RBC_velocity_.value()*t);

        // find out the sign of the gradient.
        // inlet : positive gradient
        // outlet : negative gradient
        const int sign = (z - z_inlet_ < 1e-10 ? 1 : -1);

        gradient = sign*P50_/hillExponent_.value()
            * pow(S/(1.-S),1./hillExponent_.value() - 1.)
            * beta_/sqr(1.-S);

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedGradientFvPatchField<Type>::updateCoeffs();
}

template<class Type>
void fixedPO2GradientKroghFvPatchField<Type>::write(Ostream& os) const
{
    fixedGradientFvPatchField<Type>::write(os);
    os.writeKeyword("S_inlet") << S_inlet_ << token::END_STATEMENT << nl;
    os.writeKeyword("z_inlet") << z_inlet_ << token::END_STATEMENT << nl;
    os.writeKeyword("beta") << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("P50") << P50_ << token::END_STATEMENT << nl;
    os.writeKeyword("hillExponent") << hillExponent_ << token::END_STATEMENT << nl;
    os.writeKeyword("RBC_velocity") << RBC_velocity_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
