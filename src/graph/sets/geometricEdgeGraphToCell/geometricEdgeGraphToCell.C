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

#include "geometricEdgeGraphToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(geometricEdgeGraphToCell, 0);

addToRunTimeSelectionTable(topoSetSource, geometricEdgeGraphToCell, word);

addToRunTimeSelectionTable(topoSetSource, geometricEdgeGraphToCell, istream);

}


Foam::topoSetSource::addToUsageTable Foam::geometricEdgeGraphToCell::usage_
(
    geometricEdgeGraphToCell::typeName,
    "\n    Usage: geometricEdgeGraphToCell (pt0 .. ptn)\n\n"
    "    Select the nearest cell for each of the points pt0 ..ptn\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void 
Foam::geometricEdgeGraphToCell::combine(topoSet& set, const bool add) const
{
    // cache graph points used for testing
    pointField graphPoints(graph_.nEdges()*nSamplePointsPerEdge_);

    for(int i=0; i < graph_.nEdges(); i++)
    {
        const geometricEdge& e = graph_.edge(i);
        const scalar eLength = e.length();
        
        for(int j=0; j < nSamplePointsPerEdge_; j++)
        {
            scalar s = eLength*j/(nSamplePointsPerEdge_ - 1);
            graphPoints[i*nSamplePointsPerEdge_ + j] = e.pointPosition(s);
        }
    }

    const pointField& ctrs = mesh_.cellCentres();
    forAll(ctrs, cellI)
    {
        if (isWithinMaxDistanceToGraphPoints(ctrs[cellI], graphPoints))
        {
            addOrDelete(set, cellI, add);
        }
    }
}


bool
Foam::geometricEdgeGraphToCell::isWithinMaxDistanceToGraphPoints
(
    const point& p,
    const pointField& graphPoints
) const
{
    forAll(graphPoints, i)
    {
        if (magSqr(p - graphPoints[i]) <= sqr(maxDistance_))
        {
            return true;
        }
    }

    return false;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::geometricEdgeGraphToCell::geometricEdgeGraphToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    graph_
    (
        IOobject
        (
            dict.lookup("graphDict"),
            mesh.time().constant(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    maxDistance_(readScalar(dict.lookup("maxDistance"))),
    nSamplePointsPerEdge_(dict.lookupOrDefault("samplePointsPerEdge", 100))
{}


// Construct from Istream
Foam::geometricEdgeGraphToCell::geometricEdgeGraphToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    graph_
    (
        IOobject
        (
            "graphDict",
            mesh.time().constant(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    maxDistance_(readScalar(checkIs(is))),
    nSamplePointsPerEdge_(readLabel(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::geometricEdgeGraphToCell::~geometricEdgeGraphToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::geometricEdgeGraphToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells within distance " << maxDistance_ 
            << " to the graph " << graph_.name() << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Adding cells within distance " << maxDistance_ 
            << " to the graph " << graph_.name() << endl;

        combine(set, false);
    }
}


// ************************************************************************* //

