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

#include "circularTubeGraph.H"

#include "polygonalTube.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(circularTubeGraph, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::circularTubeGraph::circularTubeGraph(const IOobject& io)
:
    geometricEdgeGraph(io),
    tubes_(nEdges())
{
    const scalarListList segmentDiameters(lookup("segmentDiameters"));

    const dictionary tubeOptions = subDict("tubeOptions");

    // sanity checks
    if (segmentDiameters.size() != nEdges())
    {
        FatalErrorIn("circularTubeGraph(const IOobject&)")
            << "The number of elements in diameters does not match the "
            << "number of edges." << nl
            << abort(FatalError);
    }

    // construct tube objects
    forAll(tubes_, eI)
    {
        tubes_.set
        (
            eI, 
            edge(eI).makeCircularTube
            (
                segmentDiameters[eI],
                tubeOptions
            )
        );
    }

    // compute ellipseAxes at the end points of each tube using the graph
    // structure
    forAll(tubes_, eI)
    {
        labelPair vv = edgeVertexIndices(eI);

        tubes_[eI].setStartAxes(vertexTubeAxes(vv.first(),  eI));
        tubes_[eI].setEndAxes  (vertexTubeAxes(vv.second(), eI));
    }

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar 
Foam::circularTubeGraph::vertexDiameter(const label vertexIndex) const
{
    scalar diameter = 0.0;

    const labelList adjacentEdges = adjacentEdgeIndices(vertexIndex);
    forAll(adjacentEdges, i)
    {
        const label eI = adjacentEdges[i];
        diameter += edgeDiameterAtVertex(vertexIndex, eI);
    }
    diameter /= adjacentEdges.size();

    return diameter;
}


Foam::scalar 
Foam::circularTubeGraph::vertexRadius(const label vertexIndex) const
{
    return 0.5*vertexDiameter(vertexIndex);
}


Foam::scalar
Foam::circularTubeGraph::edgeDiameterAtVertex
(
    const label vertexIndex,
    const label edgeIndex
) const
{
    checkVertexInEdge(vertexIndex, edgeIndex);

    labelPair vv = edgeVertexIndices(edgeIndex);
    scalar diameter = 0.0;

    if (vertexIndex == vv.first())
    {
        diameter = tubes_[edgeIndex].outerDiameter(0.0);
    }
    else if (vertexIndex == vv.second())
    {
        diameter = tubes_[edgeIndex].outerDiameter(edge(edgeIndex).length());
    }

    return diameter;
}


Foam::scalar 
Foam::circularTubeGraph::edgeRadiusAtVertex
(
    const label vertexIndex,
    const label edgeIndex
) const
{
    return 0.5*edgeDiameterAtVertex(vertexIndex, edgeIndex);
}


Foam::ellipseAxes
Foam::circularTubeGraph::vertexTubeAxes
(
    const label vertexIndex,
    const label edgeIndex
) const
{
    checkVertexInEdge(vertexIndex, edgeIndex);

    scalar minorRadius = vertexRadius(vertexIndex);
    scalar majorRadius = 
        minorRadius/Foam::sin(majorBisectionAngle(vertexIndex, edgeIndex));

    if (debug)
    {
        Info<< "For edgeIndex: " << edgeIndex
            << " and vertexIndex " << vertexIndex << ", "
            << " the ellipseAxes radii are : " 
            << minorRadius << ", " << majorRadius << nl
            << "The edge mean radius is " 
            << 0.5*tube(edgeIndex).meanOuterDiameter() << endl;
    }

    return ellipseAxes
           (
               normalToVertexPlane(vertexIndex, edgeIndex),
               majorBisectionVector(vertexIndex, edgeIndex),
               minorRadius,
               majorRadius
           );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::circularTubeGraph::~circularTubeGraph()
{}


// ************************************************************************* //
