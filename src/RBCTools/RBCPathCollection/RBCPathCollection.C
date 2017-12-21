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

#include "RBCPathCollection.H"

#include "RBCPath.H"
#include "Ostream.H"
#include "Time.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RBCPathCollection, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBCPathCollection::RBCPathCollection(const IOobject& io)
:
    IOdictionary(io),
    RBCPaths_()
{
    wordList toc = this->toc();
    forAll(toc, i)
    {
        word key = toc[i];
        if (this->isDict(key))
        {
            if (debug)
            {
                Info<< "RBCPathCollection::RBCPathCollection(const IOobject&):"
                    << " Adding RBCPath object for key " << key << "." << endl;
            }
            label index = readLabel(this->subDict(key).lookup("index"));
            RBCPaths_.insert
            (
                index, 
                autoPtr<RBCPath>
                (
                    new RBCPath(this->subDict(key))
                )
            );
            RBCIdxList_.append(index);
        }
        else
        {
            WarningIn("RBCPathCollection::RBCPathCollection(const IOobject&)")
                << "  Ignored non-dictionary key " << key << endl;
        }
    }
}


Foam::RBCPathCollection::RBCPathCollection
(
    const IOobject& io,
    const Map<autoPtr<RBCPath> >& RBCPaths    
)
:
    IOdictionary(io),
    RBCPaths_(RBCPaths)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBCPathCollection::~RBCPathCollection()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label
Foam::RBCPathCollection::size() const
{
    return RBCPaths_.size();
}


Foam::labelList
Foam::RBCPathCollection::activeIndices() const
{
    labelList activeList = labelList();
    scalar currentTime = this->time().timeOutputValue();

    forAllConstIter(Map<autoPtr<RBCPath> >, RBCPaths_, iter)
    {
        const RBCPath& cPath = *iter;
        if (cPath.active(currentTime))
        {
            activeList.append(cPath.index());
        }
    }
    return activeList;
}


Foam::labelList
Foam::RBCPathCollection::inactiveIndices() const
{
    labelList inactiveList = labelList();
    scalar currentTime = this->time().timeOutputValue();

    forAllConstIter(Map<autoPtr<RBCPath> >, RBCPaths_, iter)
    {
        const RBCPath& cPath = *iter;
        if (!cPath.active(currentTime))
        {
            inactiveList.append(cPath.index());
        }
    }
    return inactiveList;
}


Foam::labelList
Foam::RBCPathCollection::newActiveIndices() const
{
    labelList newActiveList = labelList();
    scalar currentTime = this->time().timeOutputValue();
    scalar oldTime     = this->time().timeOutputValue()
                       - this->time().deltaT().value();

    forAllConstIter(Map<autoPtr<RBCPath> >, RBCPaths_, iter)
    {
        const RBCPath& cPath = *iter;
        if (cPath.newActive(currentTime, oldTime))
        {
            newActiveList.append(cPath.index());
        }
    }

    return newActiveList;
}


Foam::labelList
Foam::RBCPathCollection::newInactiveIndices() const
{
    labelList newInactiveList = labelList();
    scalar currentTime = this->time().timeOutputValue();
    scalar oldTime     = this->time().timeOutputValue()
                       - this->time().deltaT().value();

    forAllConstIter(Map<autoPtr<RBCPath> >, RBCPaths_, iter)
    {
        const RBCPath& cPath = *iter;
        if (cPath.newInactive(currentTime, oldTime))
        {
            newInactiveList.append(cPath.index());
        }
    }

    return newInactiveList;
}


bool
Foam::RBCPathCollection::writeData(Ostream& os) const
{
    forAllConstIter(Map<autoPtr<RBCPath> >, RBCPaths_, iter)
    {
        (*iter)->write(os);
    }
    return os.good();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

const Foam::RBCPath&
Foam::RBCPathCollection::operator()(const label RBCI) const
{
    return RBCPaths_[RBCI]();
}


Foam::RBCPath&
Foam::RBCPathCollection::operator()(const label RBCI)
{
    return RBCPaths_[RBCI]();
}


const Foam::RBCPath&
Foam::RBCPathCollection::operator[](const label i) const
{
    const label RBCI = RBCIdxList_[i];
    return RBCPaths_[RBCI]();
}


Foam::RBCPath&
Foam::RBCPathCollection::operator[](const label i)
{
    const label RBCI = RBCIdxList_[i];
    return RBCPaths_[RBCI]();
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream& 
Foam::operator<<(Foam::Ostream& os, const RBCPathCollection& obj)
{
    obj.writeData(os);
    return os;
}


// ************************************************************************* //

