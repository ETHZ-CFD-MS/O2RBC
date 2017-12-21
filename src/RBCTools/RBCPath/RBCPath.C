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

#include "RBCPath.H"

#include "Ostream.H"
#include "dictionary.H"


// * * * * * * * * * * * * * Static Data Member      * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void
Foam::RBCPath::checkLists()
{
    // check that lists have the same size
    if (!((pathTimes_.size() == edges_.size()) 
       && (pathTimes_.size() == sCoords_.size())))
    {
        FatalErrorIn
        (
            "Foam::RBCPath::read()"
        )   << "The number of elements in times, edges and "
            << "sCoords for path " << index_ << " are not equal" << nl
            << exit(FatalError);
    }
    //
    // check that the lists are not empty
    if (pathTimes_.size() == 0)
    {
        FatalErrorIn
        (
            "Foam::RBCPath::read()"
        )   << "The RBCPath with index " << index_ << " is empty." 
            << exit(FatalError);
    }
    // check that the times are strictly increasing.
    for(int i=0; i < pathTimes_.size()-1; ++i)
    {
        if (pathTimes_[i] >= pathTimes_[i+1])
        {
            FatalErrorIn
            (
                "Foam::RBCPath::read()"
            )   << "In RBCPath with index " << index_ << ", the elements in "
                << "pathTimes are not strictly increasing." << nl
                << "The time " << pathTimes_[i+1] << " is not greater than " 
                << pathTimes_[i+1] << "."
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBCPath::RBCPath
(
    const dictionary& dict
)
:
    index_(readLabel(dict.lookup("index"))),
    // use SLList to read the lists from the ITStream since it supports both 
    // ascii and binary stream formats, unlike List which currently cannot 
    // read binary lists...
    pathTimes_(scalarList(SLList<scalar>(dict.lookup("times")))),
    edges_(labelList(SLList<label>(dict.lookup("edges")))),
    sCoords_(scalarList(SLList<scalar>(dict.lookup("sCoords"))))
{
    checkLists();
}


Foam::RBCPath::RBCPath
(
    const label index,
    const scalarList& pathTimes,
    const labelList& edges,
    const scalarList& sCoords
)
:
    index_(index),
    pathTimes_(pathTimes),
    edges_(edges),
    sCoords_(sCoords)
{
    checkLists();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void 
Foam::RBCPath::scaleSCoord(const label i, const scalar lambda)
{
    sCoords_[i] *= lambda;
}


void
Foam::RBCPath::setEdge(const label i, const label edgeIndex)
{
    edges_[i] = edgeIndex;
}


bool
Foam::RBCPath::active(const scalar currentTime) const
{
    // convention: the active time is defined by the interval
    //    [firstTime, lastTime)
    scalar firstTime = pathTimes_[0];
    scalar lastTime  = pathTimes_[pathTimes_.size()-1];
    if (firstTime <= currentTime && currentTime < lastTime)
    {
        return true;
    }
    return false;
}


bool
Foam::RBCPath::newActive
(
    const scalar currentTime, 
    const scalar oldTime
) const
{
    if (active(currentTime) && (!active(oldTime)))
    {
        return true;
    }
    return false;
}


bool
Foam::RBCPath::newInactive
(
    const scalar currentTime, 
    const scalar oldTime
) const
{
    if (active(oldTime) && (!active(currentTime)))
    {
        return true;
    }
    return false;
}


void
Foam::RBCPath::write(Foam::Ostream& os) const
{
    // forces to write ascii data since there are issues with binary output
    IOstream::streamFormat osFormat = os.format(IOstream::ASCII);
    os << "RBC" << index_ << nl;
    os << indent << token::BEGIN_BLOCK << incrIndent << nl;
    os.writeKeyword("index") << index_ << token::END_STATEMENT << nl;
    os.writeKeyword("times") << pathTimes_ << token::END_STATEMENT << nl;
    os.writeKeyword("edges") << edges_ << token::END_STATEMENT << nl;
    os.writeKeyword("sCoords") << sCoords_ << token::END_STATEMENT << nl;
    os << decrIndent << indent << token::END_BLOCK << endl;
    os.format(osFormat);
}


// * * * * * * * * * * * * * * IOstream Operators * * * * * * * * * * * * * * //

Foam::Ostream& 
Foam::operator<<(Ostream& os, const RBCPath& obj)
{
    obj.write(os);
    return os;
}

// ************************************************************************* //
