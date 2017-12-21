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

#include "deformableBodyGeometricState.H"

#include "dictionaryEntry.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(deformableBodyGeometricState, 0);
    defineRunTimeSelectionTable(deformableBodyGeometricState, dictionary);
    defineRunTimeSelectionTable(deformableBodyGeometricState, mesh);
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::deformableBodyGeometricState> 
Foam::deformableBodyGeometricState::New
(
    const dictionary& dict
)
{
    const word deformableBodyType = dict.lookup("type");

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(deformableBodyType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "deformableBodyGeometricState::New(const dictionary&)"
        )   << "Unknown deformable body type " << deformableBodyType
            << nl << nl
            << "Valid deformable body types : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<deformableBodyGeometricState>(cstrIter()(dict));
}


Foam::autoPtr<Foam::deformableBodyGeometricState> 
Foam::deformableBodyGeometricState::New
(
    const dictionary& dict,
    const polyMesh& mesh
)
{
    const word deformableBodyType = dict.lookup("type");

    meshConstructorTable::iterator cstrIter =
        meshConstructorTablePtr_->find(deformableBodyType);

    if (cstrIter == meshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "deformableBodyGeometricState::New(const dictionary&, "
            "const polyMesh&)"
        )   << "Unknown deformable body type " << deformableBodyType
            << nl << nl
            << "Valid deformable body types : " << endl
            << meshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<deformableBodyGeometricState>(cstrIter()(dict, mesh));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::deformableBodyGeometricState::deformableBodyGeometricState
(
    const point& center
)
:
    center_(center)
{}


Foam::deformableBodyGeometricState::deformableBodyGeometricState
(
    const dictionary& dict
)
:
    center_(dict.lookup("center"))
{}


Foam::deformableBodyGeometricState::deformableBodyGeometricState
(
    Istream& is
)
:
    center_(point::zero)
{
    read(is);
}


Foam::deformableBodyGeometricState::deformableBodyGeometricState
(
    const deformableBodyGeometricState& obj
)
:
    center_(obj.center_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Istream&
Foam::deformableBodyGeometricState::read(Istream& is)
{
    dictionary dict(is);
    center_ = dict.lookup("center");

    return is;
}


void
Foam::deformableBodyGeometricState::write(Ostream& os) const
{
    os.writeKeyword("type")   << type()  << token::END_STATEMENT << nl;
    os.writeKeyword("center") << center_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::deformableBodyGeometricState::operator=
(
    const deformableBodyGeometricState& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::deformableBodyGeometricState::operator="
                     "(const Foam::deformableBodyGeometricState&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    center_ = rhs.center_;
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

Foam::Istream&
Foam::operator>>(Istream& is, deformableBodyGeometricState& obj)
{
    obj.read(is);

    return is;
}

Foam::Ostream&
Foam::operator<<(Ostream& os, const deformableBodyGeometricState& obj)
{
    obj.write(os);

    return os;
}

// ************************************************************************* //
