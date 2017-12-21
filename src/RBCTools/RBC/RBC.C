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

#include "RBC.H"

#include "deformableBodyGeometricState.H"
#include "axisymmetricBodyGeometricState.H"

#include "dissociationCurve.H"
#include "transformField.H"


RBC::RBC
(
    const label index,
    const deformableBodyGeometricState& geometricState,
    const word& sampleMeshName,
    const scalar initialPO2
)
:
    index_(index),
    geometricStatePtr_(geometricState.clone()),
    sampleMeshName_(sampleMeshName),
    initialPO2_(initialPO2),
    inletPO2_(initialPO2)
{}


Foam::RBC::RBC
(
    const dictionary& dict
)
:
    index_(readLabel(dict.lookup("index"))),
    geometricStatePtr_
    (
        deformableBodyGeometricState::New(dict.subDict("geometricState"))
    ),
    sampleMeshName_(dict.lookup("sampleMeshName")),
    initialPO2_(readScalar(dict.lookup("initialPO2"))),
    inletPO2_(readScalar(dict.lookup("inletPO2")))
{}


Foam::RBC::RBC
(
    const RBC& obj
)
:
    index_(obj.index_),
    geometricStatePtr_(obj.geometricState().clone()),
    sampleMeshName_(obj.sampleMeshName_),
    initialPO2_(obj.initialPO2_),
    inletPO2_(obj.inletPO2_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

point
Foam::RBC::center() const
{
    return geometricStatePtr_->center();
}


void Foam::RBC::setGeometricState
(
    const deformableBodyGeometricState& geometricState
)
{
    geometricStatePtr_ = geometricState.clone();
}


void 
Foam::RBC::write(Ostream& os) const
{
    os.writeKeyword("index")      << index_      << token::END_STATEMENT << nl;
    os.writeKeyword("geometricState") << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;
    geometricState().write(os);
    os  << decrIndent << indent << token::END_BLOCK << nl;
    os.writeKeyword("sampleMeshName") << sampleMeshName_ 
                                      << token::END_STATEMENT << nl;
    os.writeKeyword("initialPO2") << initialPO2_ << token::END_STATEMENT << nl;
    os.writeKeyword("inletPO2")   << inletPO2_   << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    RBC& cRBC
)
{
    // read dictionary from Istream
    const dictionary inputDict(is);
    cRBC.index_  = readLabel(inputDict.lookup("index"));
    cRBC.setGeometricState
    (
        deformableBodyGeometricState::New(inputDict.subDict("geometricState"))()
    );
    cRBC.sampleMeshName_ = word(inputDict.lookup("sampleMeshName"));
    cRBC.initialPO2_ = readScalar(inputDict.lookup("initialPO2"));

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, RBC&)");

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const RBC& cRBC
)
{
    cRBC.write(os);

    return os;
}


// ************************************************************************* //
