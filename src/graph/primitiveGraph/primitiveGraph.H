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
    Foam::primitiveGraph

Description
    Contains the graph structure of graph data loaded by IOgraph.
    Provides basic connectivity information, such as adjacent edges and 
    vertices. 
    
    This class is templated on vertex and edge properties.

    The interface of this class uses external indices, i.e., those that are 
    read by the class IOgraph.

SourceFiles
    primitiveGraph.C

\*---------------------------------------------------------------------------*/

#ifndef primitiveGraph_H
#define primitiveGraph_H

#include "IOgraph.H"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/unordered_map.hpp>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;

// Forward declaration of friend functions and operators

template<typename VProperties, typename EProperties>
class primitiveGraph;

template<typename VProperties, typename EProperties>
Ostream& operator<<
(
    Ostream&, 
    const primitiveGraph<VProperties, EProperties>&
);

/*---------------------------------------------------------------------------*\
                         Class primitiveGraph Declaration
\*---------------------------------------------------------------------------*/

template<typename VProperties, typename EProperties>
class primitiveGraph
:
    public IOgraph
{

public:
    // typedefs 
    typedef typename boost::adjacency_list <
        boost::vecS, boost::vecS, boost::bidirectionalS, VProperties, EProperties
    > Graph;
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor VertexD;
    typedef typename GraphTraits::edge_descriptor   EdgeD;
    typedef typename GraphTraits::vertex_iterator vertex_iter;


private:
    // Private data

        //- Map of vertex index to vertex descriptor
        boost::unordered_map<label, VertexD> vertexIndexToDescriptor_;

        //- Map of edge index to edge descriptor
        boost::unordered_map<label, EdgeD> edgeIndexToDescriptor_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        primitiveGraph(const primitiveGraph&);

        //- Disallow default bitwise assignment
        void operator=(const primitiveGraph&);
        

protected:
    // Protected data

        //- Graph object
        Graph graph_;

    // Protected member functions

        //- Find vertex descriptor associated to vertex index
        VertexD findVertex(const label externalVertexIndex) const;

        //- Find edge descriptor associated to edge index
        EdgeD findEdge(const label externalEdgeIndex) const;


public:

    // Constructors

        //- Construct from IOobject
        primitiveGraph(const IOobject&);


    //- Destructor
        ~primitiveGraph();


    // Member Functions

        // Access

            //- Vertex degree
            label degree(const label externalVertexIndex) const;

            //- Return the pair of vertex descriptors for a given edge
            Pair<VertexD> edgeVertices(const label externalEdgeIndex) const;

            //- Return the vertices that are adjacent to a given vertex
            List<VertexD> adjacentVertices(const label externalVertexIndex) const;

            //- Return the edges that are adjacent to a given vertex
            List<EdgeD> adjacentEdges(const label externalVertexIndex) const;


        // Write
            bool writeData(Ostream& os) const;


    // Friend Operators

    // Ostream Operators

        friend Ostream& operator<< <VProperties, EProperties>
        (
            Ostream&, 
            const primitiveGraph<VProperties, EProperties>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "primitiveGraph.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
