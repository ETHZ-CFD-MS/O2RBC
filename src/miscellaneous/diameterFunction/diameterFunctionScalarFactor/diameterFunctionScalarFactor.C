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

#include "diameterFunctionScalarFactor.H"
#include "addToRunTimeSelectionTable.H"

#include <math.h>

namespace Foam
{
    defineTypeNameAndDebug(diameterFunctionScalarFactor, 0);

    addToRunTimeSelectionTable
    (
        diameterFunction,
        diameterFunctionScalarFactor,
        dictionary
    );
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::diameterFunctionScalarFactor::diameterFunctionScalarFactor
(
    const dictionary& dict
)
:
    diameterFunction(dict),
    factor_(readScalar(dict.lookup("factor"))),
    minDiameter_(dict.lookupOrDefault<scalar>("minDiameter", 0.0)),
    maxDiameter_(dict.lookupOrDefault<scalar>("maxDiameter", GREAT))
{}


Foam::diameterFunctionScalarFactor::diameterFunctionScalarFactor
(
    const scalar factor,
    const scalar minDiameter,
    const scalar maxDiameter
)
:
    diameterFunction(dictionary()),
    factor_(factor),
    minDiameter_(minDiameter),
    maxDiameter_(maxDiameter)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::diameterFunctionScalarFactor::operator() (const scalar diameter) const
{
    return min(max(factor_*diameter, minDiameter_), maxDiameter_);
}


void diameterFunctionScalarFactor::write(Ostream& os) const
{
    diameterFunction::write(os);
    os.writeKeyword("factor") << factor_ << token::END_STATEMENT << nl;
    os.writeKeyword("minDiameter") << minDiameter_ 
                                   << token::END_STATEMENT << nl;
    os.writeKeyword("maxDiameter") << maxDiameter_ 
                                   << token::END_STATEMENT << nl;
}


// ************************************************************************* //
