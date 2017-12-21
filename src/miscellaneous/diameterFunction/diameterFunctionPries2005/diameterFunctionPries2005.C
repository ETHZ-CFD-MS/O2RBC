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

#include "diameterFunctionPries2005.H"
#include "addToRunTimeSelectionTable.H"

#include <math.h>

namespace Foam
{
    defineTypeNameAndDebug(diameterFunctionPries2005, 0);

    addToRunTimeSelectionTable
    (
        diameterFunction,
        diameterFunctionPries2005,
        dictionary
    );
}

// * * * * * * * * * * * *  Private member functions * * * * * * * * * * * * //

Foam::scalar
Foam::diameterFunctionPries2005::ESLThickness(const scalar diameter) const
{
    // The empirical fit function is designed for lengh units of microns.
    scalar scalingFactor = 1e-6;
    scalar was = 0;
    scalar wPeak = 0;

    scalar d = diameter/scalingFactor;

    if (d <= dOff_)
    {
        was   = 0;
        wPeak = 0;
    }
    else if (d > dOff_ && d <= dCrit_)
    {
        was   = (d - dOff_)/(d + d50_ - 2*dOff_)*wMax_;
        wPeak = eAmp_*(d - dOff_)/(dCrit_ - dOff_);
    }
    else if (d > dCrit_ && d <= dTop_)
    {
        was   = (d - dOff_)/(d + d50_ - 2*dOff_)*wMax_;
        wPeak = eAmp_*std::exp(-eWidth_*(d - dCrit_));
    }
    else if (d > dTop_)
    {
        was   = (dTop_ - dOff_)/(dTop_ + d50_ - 2*dOff_)*wMax_;
        wPeak = eAmp_*std::exp(-eWidth_*(d - dCrit_));
    }

    return (was + wPeak*ePeak_) * scalingFactor;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::diameterFunctionPries2005::diameterFunctionPries2005
(
    const dictionary& dict
)
:
    diameterFunction(dict),
    dOff_  (readScalar(dict.lookup("dOff"))),
    dCrit_ (readScalar(dict.lookup("dCrit"))),
    dTop_  (readScalar(dict.lookup("dTop"))),
    d50_   (readScalar(dict.lookup("d50"))),
    eAmp_  (readScalar(dict.lookup("eAmp"))),
    eWidth_(readScalar(dict.lookup("eWidth"))),
    ePeak_ (readScalar(dict.lookup("ePeak"))),
    wMax_  (readScalar(dict.lookup("wMax")))
{}


Foam::scalar
Foam::diameterFunctionPries2005::operator() (const scalar diameter) const
{
    return (diameter - 2.*ESLThickness(diameter));
}


void diameterFunctionPries2005::write(Ostream& os) const
{
    diameterFunction::write(os);
    os.writeKeyword("dOff")   << dOff_   << token::END_STATEMENT << nl;
    os.writeKeyword("dCrit")  << dCrit_  << token::END_STATEMENT << nl;
    os.writeKeyword("dTop")   << dTop_   << token::END_STATEMENT << nl;
    os.writeKeyword("d50")    << d50_    << token::END_STATEMENT << nl;
    os.writeKeyword("eAmp")   << eAmp_   << token::END_STATEMENT << nl;
    os.writeKeyword("eWidth") << eWidth_ << token::END_STATEMENT << nl;
    os.writeKeyword("ePeak")  << ePeak_  << token::END_STATEMENT << nl;
    os.writeKeyword("wMax")   << wMax_   << token::END_STATEMENT << nl;
}


// ************************************************************************* //
