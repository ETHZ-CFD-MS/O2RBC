/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "regionDependentField.H"

#include "vascularGraphRegions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionDependentField, 1);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionDependentField::regionDependentField
(
    const IOobject& io,
    const fvMesh& mesh,
    const dictionary& dict,
    const vascularGraphRegions& regions
)
:
    field_
    (
        io,
        mesh
    ),
    regions_(regions),
    RBCValue_(),
    plasmaValue_(),
    wallValue_(),
    tissueValue_(),
    averagingMethod_(regionDependentField::ARITHMETIC)
{
    RBCValue_        = readScalar(dict.lookup("RBCValue"));
    plasmaValue_     = readScalar(dict.lookup("plasmaValue"));
    wallValue_       = readScalar(dict.lookup("wallValue"));
    tissueValue_     = readScalar(dict.lookup("tissueValue"));
    averagingMethod_ = wordToAveragingMethod(dict.lookup("averaging"));

    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionDependentField::~regionDependentField()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word
Foam::regionDependentField::averagingMethodToWord
(
    const averagingMethod& averaging
) const
{
    word averagingName;

    switch (averaging)
    {
        case regionDependentField::ARITHMETIC:
            averagingName = "arithmetic";
            break;
        case regionDependentField::HARMONIC:
            averagingName = "harmonic";
            break;
        default:
            FatalErrorIn
            (
                "Foam::regionDependentField::averagingMethodToWord(const averagingMethod&)"
            )   << "Undefined averagingMethod " << averaging << nl
                << abort(FatalError);
    }

    return averagingName;
}


regionDependentField::averagingMethod 
Foam::regionDependentField::wordToAveragingMethod(const word& averaging) const
{
    if (averaging == "arithmetic")
    {
        return regionDependentField::ARITHMETIC;
    }
    else if (averaging == "harmonic")
    {
        return regionDependentField::HARMONIC;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::regionDependentField::wordToAveragingMethod(const word&)"
        )   << "Invalid averagingMethod specifier " << averaging << nl
            << abort(FatalError);
    }

    return regionDependentField::ARITHMETIC;
}


void Foam::regionDependentField::update()
{
    // multiply by a dimensioned scalar to get the correct dimensions
    dimensionedScalar dimensionedOne("one", field_.dimensions(), 1.0);

    switch (averagingMethod_)
    {
        case regionDependentField::ARITHMETIC:
            field_ = dimensionedOne
                       * (regions_.inLumen()  * plasmaValue_
                       + regions_.inWall()   * wallValue_
                       + regions_.inTissue() * tissueValue_);
            break;
        case regionDependentField::HARMONIC:
            forAll(field_.mesh().C(), cI)
            {
                if (regions_.inTissue()[cI] > 0.0)
                {
                    if (regions_.inWall()[cI] > 0.0)
                    {
                        if (regions_.inLumen()[cI] > 0.0)
                        {
                            field_[cI] = weightedHarmonicMean
                                         (
                                            regions_.inLumen()[cI] , plasmaValue_,
                                            regions_.inWall()[cI]  , wallValue_,
                                            regions_.inTissue()[cI], tissueValue_
                                         );

                        }
                        else
                        {
                            field_[cI] = weightedHarmonicMean
                                         (
                                            regions_.inWall()[cI]  , wallValue_,
                                            regions_.inTissue()[cI], tissueValue_
                                         );
                        }
                    }
                    else
                    {
                        if (regions_.inLumen()[cI] > 0.0)
                        {
                            field_[cI] = weightedHarmonicMean
                                         (
                                            regions_.inLumen()[cI] , plasmaValue_,
                                            regions_.inTissue()[cI], tissueValue_
                                         );

                        }
                        else
                        {
                            field_[cI] = regions_.inTissue()[cI] * tissueValue_;
                        }
                    }
                }
                else
                {
                    if (regions_.inWall()[cI] > 0.0)
                    {
                        if (regions_.inLumen()[cI] > 0.0)
                        {
                            field_[cI] = weightedHarmonicMean
                                         (
                                            regions_.inLumen()[cI], plasmaValue_,
                                            regions_.inWall()[cI] , wallValue_
                                         );
                        }
                        else
                        {
                            field_[cI] = regions_.inWall()[cI] * wallValue_;
                        }
                    }
                    else
                    {
                        if (regions_.inLumen()[cI] > 0.0)
                        {
                            field_[cI] = regions_.inLumen()[cI] * plasmaValue_;
                        }
                        else
                        {
                            WarningIn
                            (
                                "Foam::regionDependentField::update()"
                            )   << "All regions have volume fraction zero or less in "
                                << "cell " << cI << "." << nl
                                << "Setting the field value to zero in this cell." << endl;

                            field_[cI] = 0.0;
                        }
                    }
                }
            }
            forAll(field_.mesh().boundary(), i)
            {
                field_.boundaryFieldRef()[i] = field_.boundaryField()[i].patchInternalField();
            }
            break;
    }
}


void Foam::regionDependentField::update(const volScalarField& inRBC)
{
    this->update();

    // for now, just use the internal field for the RBC volume fraction
    field_.primitiveFieldRef() = inRBC.primitiveField()*RBCValue_
                               + (1.0 - inRBC.primitiveField())*field_.primitiveField();
}

// writeData member function required by regIOobject
bool Foam::regionDependentField::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}

//
bool Foam::regionDependentField::write() const
{
    return field_.write();
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const regionDependentField& gf
)
{
    gf.field_.writeData(os);

    os << nl;
    os.writeKeyword("RBCValue")      << gf.RBCValue_    << token::END_STATEMENT << nl;
    os.writeKeyword("plasmaValue")   << gf.plasmaValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("wallValue")     << gf.wallValue_   << token::END_STATEMENT << nl;
    os.writeKeyword("tissueValue")   << gf.tissueValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("averingMethod") << gf.averagingMethodToWord(gf.averagingMethod_)
                                                        << token::END_STATEMENT << nl;

    // Check state of IOstream
    os.check
    (
        "Ostream& operator<<(Ostream&, "
        "const regionDependentField"
    );

    return (os);
}

// ************************************************************************* //
