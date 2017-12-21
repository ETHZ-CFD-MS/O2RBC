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

#include "RBCMover.H"

#include "dissociationCurve.H"
#include "RBCCollection.H"

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RBCMover, 0);
    defineRunTimeSelectionTable(RBCMover, meshRBCCollection);
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::RBCMover> Foam::RBCMover::New
(
    const fvMesh& mesh,
    RBCCollection& RBCCollection
)
{
    const word functionType
    (
        IOdictionary
        (
            IOobject
            (
                "RBCMoverDict",
                mesh.time().system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("type")
    );

    meshRBCCollectionConstructorTable::iterator cstrIter =
        meshRBCCollectionConstructorTablePtr_->find(functionType);

    if (cstrIter == meshRBCCollectionConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "RBCMover::New(const fvMesh&, const RBCCollection&)"
        )   << "Unknown RBCMover type " << functionType
            << nl << nl
            << "Valid RBCMover types : " << endl
            << meshRBCCollectionConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<RBCMover>(cstrIter()(mesh, RBCCollection));
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

scalar
RBCMover::computePO2FromFieldNameAndValue
(
    const scalar value,
    const word& fieldName
) const
{
    scalar PO2Value = 0.0;
    if (fieldName == "PO2")
    {
        PO2Value = value;
    }
    else if (fieldName == "Hb")
    {
        PO2Value = RBCCollection_.getDissociationCurve().equilibriumPO2(value);
    }
    else
    {
        FatalErrorIn
        (
            "Foam::RBCMover::computePO2FromFieldNameAndValue(const RBCPath&) const:"
        ) << "Invalid field name " << fieldName << ". " << nl
          << "Supported field names for the RBC insertion are PO2 and Hb."
          << abort(FatalError);
    }
    return PO2Value;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

RBCMover::RBCMover
(
    const fvMesh& mesh,
    RBCCollection& RBCCollection
)
:
    IOdictionary
    (
        IOobject
        (
            "RBCMoverDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    RBCCollection_(RBCCollection)
{
    // intentionally blank
}

RBCMover::~RBCMover()
{
    // intentionally blank
}


bool 
RBCMover::writeData
(
    Ostream& os
) const
{
    return os.good();
}


