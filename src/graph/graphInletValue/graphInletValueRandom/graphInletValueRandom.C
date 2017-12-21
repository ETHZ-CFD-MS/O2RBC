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

#include "graphInletValueRandom.H"
#include "addToRunTimeSelectionTable.H"

#include "geometricEdgeGraph.H"

namespace Foam
{
    defineTypeNameAndDebug(graphInletValueRandom, 0);

    addToRunTimeSelectionTable
    (
        graphInletValue,
        graphInletValueRandom,
        dictionary
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::graphInletValueRandom::graphInletValueRandom
(
    const dictionary& dict,
    const geometricEdgeGraph& graph
)
:
    graphInletValue(dict, graph),
    inletRandomVars_()
{
    if (dict.found("uniformValue"))
    {
        randomVariable randomVar(dict.subDict("uniformValue"));
        inletRandomVars_ = Map<randomVariable>(graph_.nEdges());
        forAll(graph_.edgeIndices(), eI)
        {
            inletRandomVars_.set(eI, randomVar);
        }
    }
    else
    {
        inletRandomVars_ = Map<randomVariable>(dict.lookup("values"));
    }
    forAll(inletRandomVars_.toc(), eI)
    {
        assignedEdges_.append(eI);
    }
    checkLeafEdges();
}


Foam::scalar
Foam::graphInletValueRandom::inletValue
(
    const label edgeIndex,
    const scalar time
) const
{
    if (!inletRandomVars_.found(edgeIndex))
    {
        FatalErrorIn
        (
            "Foam::graphInletValueRandom::inletValue(const label, const scalar)"
        ) << "No inlet value specified for edge " << edgeIndex << "." << endl
          << abort(FatalError);
    }
    return inletRandomVars_[edgeIndex].generate();
}


void graphInletValueRandom::write(Ostream& os) const
{
    graphInletValue::write(os);
    os.writeKeyword("values") << inletRandomVars_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Istream&
Foam::operator>>(Istream& is, graphInletValueRandom& obj)
{
    dictionary dict(is);
    dict.lookup("values") >> obj.inletRandomVars_;

    return is;
}


Foam::Ostream&
Foam::operator<<(Ostream& os, const graphInletValueRandom& obj)
{
    obj.write(os);

    return os;
}


// ************************************************************************* //
