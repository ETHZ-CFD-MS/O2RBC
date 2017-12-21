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

#include "graphInletValueConstant.H"
#include "addToRunTimeSelectionTable.H"

#include "geometricEdgeGraph.H"

namespace Foam
{
    defineTypeNameAndDebug(graphInletValueConstant, 0);

    addToRunTimeSelectionTable
    (
        graphInletValue,
        graphInletValueConstant,
        dictionary
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::graphInletValueConstant::graphInletValueConstant
(
    const dictionary& dict,
    const geometricEdgeGraph& graph
)
:
    graphInletValue(dict, graph),
    inletValues_()
{
    if (dict.found("uniformValue"))
    {
        scalar uniformValue = readScalar(dict.lookup("uniformValue"));
        inletValues_ = Map<scalar>(graph_.nEdges());
        forAll(graph_.edgeIndices(), eI)
        {
            inletValues_.set(eI, uniformValue);
        }
    }
    else
    {
        inletValues_ = Map<scalar>(dict.lookup("values"));
    }
    forAll(inletValues_.toc(), eI)
    {
        assignedEdges_.append(eI);
    }
    checkLeafEdges();
}


Foam::scalar
Foam::graphInletValueConstant::inletValue
(
    const label edgeIndex,
    const scalar time
) const
{
    if (!inletValues_.found(edgeIndex))
    {
        FatalErrorIn
        (
            "Foam::graphInletValueConstant::inletValue(const label, const scalar)"
        ) << "No inlet value specified for edge " << edgeIndex << "." << endl
          << abort(FatalError);
    }
    return inletValues_[edgeIndex];
}


void graphInletValueConstant::write(Ostream& os) const
{
    graphInletValue::write(os);
    os.writeKeyword("values") << inletValues_ << token::END_STATEMENT << nl;
}



// ************************************************************************* //
