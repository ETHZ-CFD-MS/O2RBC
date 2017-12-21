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
    Foam::IOgraph

Description
    Provides IO functionality for a graph defined using an adjacency list.

    This class also supports vertex indices and edge indices that are
    arbitrary integers.

SourceFiles
    IOgraph.C

\*---------------------------------------------------------------------------*/

#ifndef IOgraph_H
#define IOgraph_H

#include "IOdictionary.H"
#include "IOList.H"
#include "labelIOList.H"
#include "labelListIOList.H"
#include "labelList.H"
#include "labelPair.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;

// Forward declaration of friend functions and operators
class IOgraph;
Ostream& operator<<(Ostream&, const IOgraph&);

/*---------------------------------------------------------------------------*\
                         Class IOgraph Declaration
\*---------------------------------------------------------------------------*/

class IOgraph
:
    public IOdictionary
{
    // Private data

        //- List of external vertex indices
        labelList vertexIndices_;

        //- List of external edge indices
        labelList edgeIndices_;

        //- Graph adjacency list
        labelPairList adjacencyList_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        IOgraph(const IOgraph&);

        //- Disallow default bitwise assignment
        void operator=(const IOgraph&);


public:

    // Constructors

        //- Construct from IOobject
        IOgraph(const IOobject&);


    //- Destructor
    ~IOgraph();


    // Member Functions

        // Database

            //- Return the local graph directory (constant/meshSubDir)
            fileName graphDir() const;

        // Access
            const labelList& vertexIndices() const
            {
                return vertexIndices_;
            }

            const labelList& externalVertexIndices() const
            {
                return vertexIndices_;
            }

            const labelList& edgeIndices() const
            {
                return edgeIndices_;
            }

            const labelList& externalEdgeIndices() const
            {
                return edgeIndices_;
            }

            const labelPairList& adjacencyList() const
            {
                return adjacencyList_;
            }

            //- Return the number of vertices
            label nVertices() const;

            //- Return the number of edges
            label nEdges() const;

            
        // Write
            
            //- Write the underlying dictionary
            void write(Ostream& os) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const IOgraph&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
