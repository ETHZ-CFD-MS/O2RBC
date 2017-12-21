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
    Foam::circularTubeGraph

Description
    Abstract class for graphs in three-dimensional space that have a circular 
    tube around their edges.

    The diameters correspond to the outer tube diameters.

SourceFiles
    circularTubeGraph.C

\*---------------------------------------------------------------------------*/

#ifndef circularTubeGraph_H
#define circularTubeGraph_H

#include "geometricEdgeGraph.H"

#include "circularTube.H"
#include "ellipseAxes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class circularTubeGraph Declaration
\*---------------------------------------------------------------------------*/

class circularTubeGraph
:
    public geometricEdgeGraph
{
    // Private data
    
        //- List of tubes
        PtrList<circularTube> tubes_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        circularTubeGraph(const circularTubeGraph&);

        //- Disallow default bitwise assignment
        void operator=(const circularTubeGraph&);

public:

    // Runtime type information
    ClassName("circularTubeGraph");

    // Constructors

        //- Construct from IOobject
        circularTubeGraph(const IOobject&);


    //- Destructor
        ~circularTubeGraph();


    // Member Functions

        // Access

            //- Constant access to tube elements
            const circularTube& tube(const label edgeIndex) const
            {
                return tubes_[edgeIndex];
            }

        // Geometric queries

            //- Return the averaged tube diameter at a given vertex.
            scalar vertexDiameter(const label vertexIndex) const;

            //- Return the averaged tube radius at a given vertex.
            scalar vertexRadius(const label vertexIndex) const;

            //- Return the tube diameter of an edge at a given vertex
            scalar edgeDiameterAtVertex
            (
                const label vertexIndex, 
                const label edgeIndex
            ) const;

            //- Return the tube radius of an edge at a given vertex
            scalar edgeRadiusAtVertex
            (
                const label vertexIndex, 
                const label edgeIndex
            ) const;

            //- Return the tube axes at a given vertex.
            //  The axes are guaranteed to be normal to the given edge.
            ellipseAxes vertexTubeAxes
            (
                const label vertexIndex,
                const label edgeIndex
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
