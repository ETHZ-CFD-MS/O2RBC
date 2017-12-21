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

#include "error.H"
#include "graphInletValue.H"

#include "geometricEdgeGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(graphInletValue, 0);
    defineRunTimeSelectionTable(graphInletValue, dictionary);
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::graphInletValue> Foam::graphInletValue::New
(
    const dictionary& dict,
    const geometricEdgeGraph& graph

)
{
    const word functionType = dict.lookup("type");

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(functionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "graphInletValue::New(const word&, "
            "const dictionary&)"
        )   << "Unknown graph inlet value type " << functionType
            << nl << nl
            << "Valid graph inlet value types : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<graphInletValue>(cstrIter()(dict, graph));
}


// * * * * * * * * * * *  Protected member functions * * * * * * * * * * * * //

void
Foam::graphInletValue::checkLeafEdges() const
{
    forAll(graph_.vertexIndices(), vI)
    {
        if (graph_.vertexDegree(vI) == 1)
        {
            labelList edges = graph_.adjacentEdgeIndices(vI);
            label eI = edges[0];
            if (findIndex(assignedEdges_, eI) < 0)
            {
                WarningIn
                (
                    "Foam::graphInletValue::checkLeafEdges() const"
                ) << "The edge " << eI << " has no inlet value specified "
                  << "although it is a leaf edge." << endl;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::graphInletValue::graphInletValue
(
    const dictionary& dict,
    const geometricEdgeGraph& graph
)
:   graph_(graph),
    fieldName_(dict.lookup("field")),
    assignedEdges_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void graphInletValue::write(Ostream& os) const
{
    os.writeKeyword("type")   << type()     << token::END_STATEMENT << nl;
    os.writeKeyword("field")  << fieldName_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream&
Foam::operator<<(Ostream& os, const graphInletValue& obj)
{
    obj.write(os);

    return os;
}


// ************************************************************************* //
