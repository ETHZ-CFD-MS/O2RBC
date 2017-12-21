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
    Foam::tissueGraphRegions

Description
    Class for the computation of the tissue regions closest to a tube graph.

SourceFiles
    tissueGraphRegions.C

\*---------------------------------------------------------------------------*/

#ifndef tissueGraphRegions_H
#define tissueGraphRegions_H

#include "fvCFD.H"

#include "labelIOField.H"

namespace Foam
{

// Forward declaration of classes
class Ostream;
class circularTubeGraph;

// Forward declaration of friend functions and operators
class tissueGraphRegions;
Ostream& operator<<(Ostream&, const tissueGraphRegions&);


/*---------------------------------------------------------------------------*\
                         Class tissueGraphRegions Declaration
\*---------------------------------------------------------------------------*/

class tissueGraphRegions
{
    // Private data

        //- Reference to Eulerian mesh
        const fvMesh& mesh_;

        //- Reference to graph
        const circularTubeGraph& graph_;

        //- Field with the indices of the nearest edge
        labelIOField nearestEdgeIndices_;

        //- Field with the scalar indices of the nearest edge
        //  Used for visualization.
        volScalarField nearestEdgeScalarIndices_;

        //- number of sample points for each graph edge
        label nSamplePointsPerEdge_;

        //- Index mapping for edges that cross the domain bounding box
        //  multiple times.
        //  A key-value pair of the map has the structure
        //  eI, [(s1, e1New), (s2, e2New), ...]
        //  This means that all points with s < s1 keep the index eI.
        //  All points with s1 < s < s2 get the index e1New etc.
        Map<List<Tuple2<scalar, label> > > edgeSegmentIndexMap_;

        
    // Private Member Functions
    
        //- Return whether the given edge overlaps the Eulerian mesh.
        //  Takes a local index as argument.
        bool edgeOverlapsMesh(const label edgeIndex) const;

        //- List of edges that overlap the Eulerian mesh.
        //  Returns local indices.
        const labelList overlappingEdgeIndices() const;

        //- Return the index of the nearest edge to the given point
        label nearestEdgeIndex
        (
            const point& p,
            const DynamicList<point>& graphPoints,
            const DynamicList<label>& graphPointsEdgeIndices
        ) const;

        //- Return the maximal edge segment index
        label maximumEdgeSegmentIndex() const;

        //- Disallow default bitwise copy construct
        tissueGraphRegions(const tissueGraphRegions&);

        //- Disallow default bitwise assignment
        void operator=(const tissueGraphRegions&);


public:

    // Runtime type information
    ClassName("tissueGraphRegions");

    // Constructors

        //- Construct by reading the vascular regions and the associated edge
        //  meshes
        tissueGraphRegions
        (
            const fvMesh& mesh, 
            const circularTubeGraph& 
        );


    //- Destructor
    ~tissueGraphRegions();


    // Member Functions

        //- Calculate regions
        void calculate();

        //- Constant access to the indices of the nearest edges
        const labelField& nearestEdgeIndices() const
        {
            return nearestEdgeIndices_;
        }

        //- Sum a scalar field on each tissue region
        Map<scalar> sumOnTissueRegions(const scalarField&) const;
        
        //- Sum a scalar field on each tissue region
        Map<scalar> sumOnTissueRegions(const volScalarField&) const;

        //- Average a scalar field on each tissue region with given integration
        //  weights.
        Map<scalar> averageOnTissueRegions
        (
            const volScalarField& f,
            const scalarField& w
        ) const;

        //- Average a scalar field on each tissue region with cell volumes as
        //  integration weights.
        Map<scalar> averageOnTissueRegions
        (
            const volScalarField& f
        ) const;

        //- Write the fields computed by this class
        void writeFields() const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const tissueGraphRegions&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
