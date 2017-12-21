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

Class
    Foam::geometricEdgeGraphToCell

Description
    A topoSetSource to select cells that are within a given distance to
    vertices on a geometricGraphObject.

SourceFiles
    geometricEdgeGraphToCell.C

\*---------------------------------------------------------------------------*/

#ifndef geometricEdgeGraphToCell_H
#define geometricEdgeGraphToCell_H

#include "topoSetSource.H"

#include "geometricEdgeGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class geometricEdgeGraphToCell Declaration
\*---------------------------------------------------------------------------*/

class geometricEdgeGraphToCell
:
    public topoSetSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- graph object
        const geometricEdgeGraph graph_;

        //- maximal distance from the graph for cell selection
        scalar maxDistance_;

        //- number of sample points for each graph edge
        label nSamplePointsPerEdge_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;

        bool isWithinMaxDistanceToGraphPoints
        (
            const point&,
            const pointField&
        ) const;


public:

    //- Runtime type information
    TypeName("geometricEdgeGraphToCell");

    // Constructors

        //- Construct from dictionary
        geometricEdgeGraphToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        geometricEdgeGraphToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    //- Destructor
    virtual ~geometricEdgeGraphToCell();


    // Member Functions

        virtual sourceType setType() const
        {
            return CELLSETSOURCE;
        }

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

