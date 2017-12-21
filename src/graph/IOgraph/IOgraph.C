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

#include "IOgraph.H"

#include "Time.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOgraph::IOgraph(const IOobject& io)
:
    IOdictionary(io),
    vertexIndices_(lookup("vertexIndices")),
    edgeIndices_(lookup("edgeIndices")),
    adjacencyList_(lookup("adjacencyList"))
{
    // sanity checks
    if (edgeIndices_.size() != adjacencyList_.size())
    {
        FatalErrorIn
        (
            "IOgraph::IOgraph(const IOobject& io)"
        )   << "  edgeIndices and adjacencyList do not have the same size" << nl
            << abort(FatalError);
    }
    labelList duplicates;
    duplicateOrder(vertexIndices_, duplicates);
    if (duplicates.size() > 0)
    {
        FatalErrorIn
        (
            "IOgraph::IOgraph(const IOobject& io)"
        )   << "  There are duplicates in vertexIndices" << nl
            << abort(FatalError);
    }
    duplicateOrder(edgeIndices_, duplicates);
    if (duplicates.size() > 0)
    {
        FatalErrorIn
        (
            "IOgraph::IOgraph(const IOobject& io)"
        )   << "  There are duplicates in edgeIndices" << nl
            << abort(FatalError);
    }
    forAll(adjacencyList_, i)
    {
        labelPair pairI = adjacencyList_[i];
        bool found1 = findIndex(vertexIndices_, pairI.first()) >= 0 ? true : false;
        bool found2 = findIndex(vertexIndices_, pairI.second()) >= 0 ? true : false;
        if (!(found1 && found2))
        {
            FatalErrorIn
            (
                "IOgraph::IOgraph(const IOobject& io)"
            )   << "  An index in the adjacency list was not found in the "
                << " vertex indices." << nl
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IOgraph::~IOgraph()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label
Foam::IOgraph::nVertices() const
{
    return vertexIndices_.size();
}


Foam::label
Foam::IOgraph::nEdges() const
{
    return edgeIndices_.size();
}


void
Foam::IOgraph::write(Ostream& os) const
{
    dictionary::write(os);
}



// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::IOgraph::operator=(const IOgraph& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::IOgraph::operator=(const Foam::IOgraph&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //
Foam::Ostream&
Foam::operator<<(Ostream& os, const IOgraph& graph)
{
    graph.write(os);
    return os;
}


// ************************************************************************* //
