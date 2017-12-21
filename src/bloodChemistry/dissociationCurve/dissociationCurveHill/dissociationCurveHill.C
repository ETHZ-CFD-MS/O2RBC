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

#include "dissociationCurveHill.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(dissociationCurveHill, 0);

    addToRunTimeSelectionTable
    (
        dissociationCurve,
        dissociationCurveHill,
        dictionary
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dissociationCurveHill::dissociationCurveHill
(
    const dictionary& dict
)
:
    dissociationCurve(dict),
    P50_(readScalar(dict.lookup("P50"))),
    exponent_(readScalar(dict.lookup("exponent")))
{}


Foam::dissociationCurveHill::dissociationCurveHill
(
    const scalar P50,
    const scalar exponent
)
:
    dissociationCurve(),
    P50_(P50),
    exponent_(exponent)
{}


Foam::dissociationCurveHill::dissociationCurveHill
(
    const dissociationCurveHill& dc
)
:
    dissociationCurve(dc),
    P50_(dc.P50_),
    exponent_(dc.exponent_)
{}


Foam::scalar 
Foam::dissociationCurveHill::equilibriumPO2
(
    const scalar& S
) const
{
    return P50_*pow(S/(1. - S),1./exponent_);
}


Foam::dimensionedScalar 
Foam::dissociationCurveHill::equilibriumPO2
(
    const dimensionedScalar& S
) const
{
    return P50_*pow(S/(1. - S),1./exponent_);
}


Foam::tmp<Foam::scalarField> 
Foam::dissociationCurveHill::equilibriumPO2
(
    const scalarField& S
) const
{
    tmp<scalarField> tP
    (
        new scalarField(P50_*pow(S/(1. - S),1./exponent_))
    );
    return tP;
}


Foam::tmp<Foam::volScalarField> 
Foam::dissociationCurveHill::equilibriumPO2
(
    const volScalarField& S
) const
{
    tmp<volScalarField> tP
    (
        new volScalarField(P50_*pow(S/(1. - S),1./exponent_))
    );
    return tP;
}


Foam::scalar 
Foam::dissociationCurveHill::equilibriumS
(
    const scalar& P
) const
{
    return pow(P,exponent_)/(pow(P,exponent_) + pow(P50_,exponent_));
}


Foam::dimensionedScalar 
Foam::dissociationCurveHill::equilibriumS
(
    const dimensionedScalar& P
) const
{
    return pow(P,exponent_)/(pow(P,exponent_) + pow(P50_,exponent_));
}


Foam::tmp<Foam::scalarField> 
Foam::dissociationCurveHill::equilibriumS
(
    const scalarField& P
) const
{
    tmp<scalarField> tS
    (
        new scalarField(pow(P,exponent_)/(pow(P,exponent_) + pow(P50_,exponent_)))
    );
    return tS;
}


Foam::tmp<Foam::volScalarField> 
Foam::dissociationCurveHill::equilibriumS
(
    const volScalarField& P
) const
{
    tmp<volScalarField> tS
    (
        new volScalarField(pow(P,exponent_)/(pow(P,exponent_) + pow(P50_,exponent_)))
    );
    return tS;
}


void Foam::dissociationCurveHill::write(Ostream& os) const
{
    dissociationCurve::write(os);
    os.writeKeyword("P50")     << P50_     << token::END_STATEMENT << nl;
    os.writeKeyword("exponent") << exponent_ << token::END_STATEMENT << nl;
}


Foam::tmp<Foam::volScalarField> 
Foam::dissociationCurveHill::reactionTerm
(
    const volScalarField& PO2,
    const volScalarField& S,
    const volScalarField& inRBC
) const
{
    tmp<volScalarField> tf
    (
        new volScalarField
        (  
            S - (inRBC - S)*pow(PO2/P50_, exponent_)
        )
    );
    return tf;
}


Foam::tmp<Foam::volScalarField> 
Foam::dissociationCurveHill::reactionTerm
(
    const volScalarField& PO2,
    const volScalarField& S
) const
{
    tmp<volScalarField> tf
    (
        new volScalarField
        (  
            S - (1.0 - S)*pow(PO2/P50_, exponent_)
        )
    );
    return tf;
}


Foam::tmp<Foam::volScalarField> 
Foam::dissociationCurveHill::ddP
(
    const volScalarField& PO2,
    const volScalarField& S,
    const volScalarField& inRBC
) const
{
    tmp<volScalarField> tddP
    (
        new volScalarField
        (  
            (S - inRBC)*exponent_ * pow(PO2, exponent_ - 1.0)/pow(P50_, exponent_)
        )
    );
    return tddP;
}


Foam::tmp<Foam::volScalarField> 
Foam::dissociationCurveHill::ddS
(
    const volScalarField& PO2,
    const volScalarField& S
) const
{
    tmp<volScalarField> tddS
    (
        new volScalarField
        (
            1.0 + pow(PO2/P50_, exponent_)
        )
    );
    return tddS;
}


// * * * * * * * * * * * * * * IOstream Operators * * * * * * * * * * * * * * //

Foam::Ostream&
Foam::operator<<(Ostream& os, const dissociationCurveHill& obj)
{
    obj.write(os);

    os.check("Ostream& operator<<(Ostream&, "
             "const dissociationCurveHill&");

    return os;
}


// ************************************************************************* //
