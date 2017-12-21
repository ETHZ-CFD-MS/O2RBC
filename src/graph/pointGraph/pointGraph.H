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
    Foam::pointGraph

Description
    Graph structure with vertices that correspond to 3D points.

    This class stores no geometric information about edges.

    The interface of this classes uses internal indices, i.e., the indices
    span the range between zero and the number of vertex/edges minus one.

SourceFiles
    pointGraph.C

\*---------------------------------------------------------------------------*/

#ifndef pointGraph_H
#define pointGraph_H

#include "primitiveGraph.H"

#include "point.H"
#include "pointField.H"
#include "scalarList.H"
#include "vectorList.H"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/unordered_map.hpp>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //- Vertex properties for pointGraph
    struct pointGraphVProp
    {
        label internalIndex;
        label externalIndex;
        point position;
    };

    //- Edge properties for pointGraph
    struct pointGraphEProp
    {
        label internalIndex;
        label externalIndex;
    };

// Forward declaration of classes
class Ostream;

// Forward declaration of friend functions and operators
class pointGraph;
Ostream& operator<<(Ostream&, const pointGraph&);

/*---------------------------------------------------------------------------*\
                         Class pointGraph Declaration
\*---------------------------------------------------------------------------*/

class pointGraph
:
    public primitiveGraph<pointGraphVProp,pointGraphEProp>
{
    // Private data

        //- Vertex positions
        vectorList vertexPositions_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointGraph(const pointGraph&);

        //- Disallow default bitwise assignment
        void operator=(const pointGraph&);


public:

    // Runtime type information
    ClassName("pointGraph");

    // Constructors

        //- Construct from IOobject
        pointGraph(const IOobject&);


    //- Destructor
        ~pointGraph();


    // Member Functions

        // Access

            //- Return vertex position
            const point& vertexPosition(const label vertexIndex) const;
        
            //- Return indices of adjacent vertices to an edge
            labelPair edgeVertexIndices(const label edgeIndex) const;

            //- Return indices of adjacent vertices to a vertex    
            labelList adjacentVertexIndices(const label vertexIndex) const;

            //- Return indices of adjacent edges to a vertex    
            labelList adjacentEdgeIndices(const label vertexIndex) const;

            //- Return positions of adjacent vertices to an edge
            Pair<point> edgeVertexPositions(const label edgeIndex) const;

            //- For a given vertex, return an adjacent edge index which is different
            //  from the provided edge index. If there is no different edge, it
            //  returns -1.
            label otherAdjacentEdgeIndex
            (
                const label vertexIndex, 
                const label edgeIndex
            ) const;

            //- For a given edge, return the vertex index which is different
            //  from the given vertex index.
            label otherEdgeVertexIndex
            (
                const label vertexIndex, 
                const label edgeIndex
            ) const;

            label vertexDegree(const label vertexIndex) const;

            //- Whether the given edge is a leaf
            bool isLeafEdge(const label edgeIndex) const;

            label externalToInternalVertexIndex
            (
                const label externalVertexIndex
            ) const;

            label externalToInternalEdgeIndex
            (
                const label externalEdgeIndex
            ) const;

        // Check
        
            //- Check if the vertex belongs to the given edge.
            void checkVertexInEdge
            (
                const label vertexIndex, 
                const label edgeIndex
            ) const;

        // Write
            bool writeData(Ostream& os) const;


    // Ostream Operators

        friend Ostream& operator<<(Ostream&, const pointGraph&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
