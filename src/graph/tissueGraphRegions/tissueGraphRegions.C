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

#include "tissueGraphRegions.H"

#include "circularTubeGraph.H"
#include "meshToMeshMoving.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tissueGraphRegions, 0);

    //- combination operator for scalar sums on tissue graph regions
    class sumDefinedOp
    {
    public:

        void operator()(scalar& x, const scalar& y) const
        {
            const scalar unsetVal(-VGREAT);

            if (x == unsetVal && y != unsetVal)
            {
                x = y;
            }
            else if (x != unsetVal && y != unsetVal)
            {
                x += y;
            }
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool
Foam::tissueGraphRegions::edgeOverlapsMesh(const label edgeIndex) const
{
    pointField tubeSkeleton = graph_.tube(edgeIndex).skeletonPoints();
    const boundBox& meshBb = mesh_.bounds();

    forAll(tubeSkeleton, i)
    {
        if (meshBb.contains(tubeSkeleton[i]))
        {
            return true;
        }
    }

    return false;
}


const Foam::labelList
Foam::tissueGraphRegions::overlappingEdgeIndices() const
{
    labelList indexList;

    forAll(graph_.edgeIndices(), eI)
    {
        if (edgeOverlapsMesh(eI))
        {
            indexList.append(eI);
        }
    }

    return indexList;
}


label Foam::tissueGraphRegions::nearestEdgeIndex
(
    const point& p,
    const DynamicList<point>& graphPoints,
    const DynamicList<label>& graphPointsEdgeIndices
) const
{
    scalar minDistanceSqr(VGREAT);
    label nearestEdgeIndex(-1);
    forAll(graphPoints, i)
    {
        scalar distanceSqr = magSqr(p - graphPoints[i]);
        if (distanceSqr < minDistanceSqr)
        {
            minDistanceSqr = distanceSqr;
            nearestEdgeIndex = graphPointsEdgeIndices[i];
        }
    }
    return nearestEdgeIndex;
}


label Foam::tissueGraphRegions::maximumEdgeSegmentIndex() const
{
    label maxSegmentIndex = graph_.nEdges();
    forAll(edgeSegmentIndexMap_.toc(), i)
    {
        label eI = edgeSegmentIndexMap_.toc()[i];
        List<Tuple2<scalar, label> > segmentMap = edgeSegmentIndexMap_[eI];
        forAll(segmentMap, sI)
        {
            Tuple2<scalar, label> tup = segmentMap[sI];
            maxSegmentIndex = max(tup.second(), maxSegmentIndex);
        }
    }
    return maxSegmentIndex;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tissueGraphRegions::tissueGraphRegions
(
    const fvMesh& mesh, 
    const circularTubeGraph& graph
)
:
    mesh_(mesh),
    graph_(graph),
    nearestEdgeIndices_
    (
        IOobject
        (
            "nearestEdgeIndices",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_.nCells()
    ),
    nearestEdgeScalarIndices_
    (
        IOobject
        (
            "nearestEdgeScalarIndices",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimless
    ),
    nSamplePointsPerEdge_(100),
    edgeSegmentIndexMap_()
{
    IOdictionary dict
    (
        IOobject
        (
            "tissueGraphRegionsDict",
            mesh_.time().caseSystem(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dictionary()
    );
    if (dict.found("edgeSegmentIndexMap"))
    {
        edgeSegmentIndexMap_ = Map<List<Tuple2<scalar, label> > >
                               (
                                    dict.lookup("edgeSegmentIndexMap")
                               );
    }
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tissueGraphRegions::~tissueGraphRegions()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
Foam::tissueGraphRegions::calculate()
{
    Info<< "Computing the indices of the nearest edges" << endl;

    // cache graph points used for testing
    DynamicList<point> graphPoints;
    DynamicList<label> graphPointsEdgeIndices;

    for(int eI=0; eI < graph_.nEdges(); eI++)
    {
        const geometricEdge& e = graph_.edge(eI);
        const scalar eLength = e.length();
        
        for(int j=0; j < nSamplePointsPerEdge_; j++)
        {
            scalar s = eLength*j/(nSamplePointsPerEdge_ - 1);
            point p = e.pointPosition(s);
            if (mesh_.bounds().contains(p))
            {
                graphPoints.append(p);
                // uses edge segment index mapping if a corresponding key is
                // present
                label segmentI = eI;
                if (edgeSegmentIndexMap_.found(eI))
                {
                    List<Tuple2<scalar, label> > segmentMap = edgeSegmentIndexMap_[eI];
                    forAll(segmentMap, sI)
                    {
                        Tuple2<scalar, label> tup = segmentMap[sI];
                        if (s > tup.first())
                        {
                            segmentI = tup.second();
                        }
                    }
                }
                graphPointsEdgeIndices.append(segmentI);
            }
        }
    }

    const pointField& ctrs = mesh_.cellCentres();
    forAll(ctrs, cellI)
    {
        nearestEdgeIndices_[cellI] = 
            nearestEdgeIndex
            (
                ctrs[cellI], 
                graphPoints, 
                graphPointsEdgeIndices
            );
        nearestEdgeScalarIndices_.primitiveFieldRef()[cellI] = nearestEdgeIndices_[cellI];
    }
    // set values on patch to be equal to the neighboring internal field.
    forAll(mesh_.boundary(), i)
    {
        nearestEdgeScalarIndices_.boundaryFieldRef()[i]
            = nearestEdgeScalarIndices_.boundaryField()[i].patchInternalField();
    }
}


Foam::Map<scalar> 
Foam::tissueGraphRegions::sumOnTissueRegions
(
    const scalarField& f
) const
{
    if (f.size() != mesh_.nCells())
    {
        FatalErrorIn("Foam::tissueGraphRegions::sumOnTissueRegions"
                     "(const volScalarField&)")
                << "the size of given field (" << f.size() << ")"
                << " does not match the number of cells (" << mesh_.nCells() << ")"
                << abort(FatalError);
    }
    // first do the summation with a list for easier parallel functionality
    const scalar unsetVal(-VGREAT);
    scalarList fSums(maximumEdgeSegmentIndex() + 1, unsetVal);
    forAll(mesh_.cells(), cellI)
    {
        label segmentI = nearestEdgeIndices_[cellI];
        if (fSums[segmentI] == unsetVal)
        {
            fSums[segmentI] = 0.0;
        }
        fSums[segmentI] += f[cellI];
    }

    Pstream::listCombineGather(fSums, sumDefinedOp());
    Pstream::listCombineScatter(fSums);
    // from the list values, create a map with all set values of the sum
    Map<scalar> fSumsMap;
    forAll(fSums, segmentI)
    {
        if (fSums[segmentI] != unsetVal)
        {
            fSumsMap.insert(segmentI, fSums[segmentI]);
        }
    }
    return fSumsMap;
}


Foam::Map<scalar> 
Foam::tissueGraphRegions::sumOnTissueRegions
(
    const volScalarField& f
) const
{
    if (f.mesh() != mesh_)
    {
        FatalErrorIn("Foam::tissueGraphRegions::sumOnTissueRegions"
                     "(const volScalarField&)")
                << "different mesh for field " << f.name()
                << abort(FatalError);
    }
    return sumOnTissueRegions(f.internalField());
}


Foam::Map<scalar> 
Foam::tissueGraphRegions::averageOnTissueRegions
(
    const volScalarField& f,
    const scalarField& w
) const
{
    if (w.size() != mesh_.nCells())
    {
        FatalErrorIn("Foam::tissueGraphRegions::sumOnTissueRegions"
                     "(const volScalarField&)")
                << "the size of given weigths (" << w.size() << ")"
                << " does not match the number of cells (" << mesh_.nCells() << ")"
                << abort(FatalError);
    }
    Map<scalar> fieldSum = sumOnTissueRegions(f.internalField()*w);
    Map<scalar> weightSum = sumOnTissueRegions(w);
    Map<scalar> fAverage;
    forAll(fieldSum.toc(), i)
    {
        label edgeI = fieldSum.toc()[i];
        if (weightSum[edgeI] != 0.0)
        {
            fAverage.insert(edgeI, fieldSum[edgeI]/weightSum[edgeI]);
        }
        else
        {
            fAverage.insert(edgeI, 0);
        }
    }
    return fAverage;
}


Foam::Map<scalar> 
Foam::tissueGraphRegions::averageOnTissueRegions
(
    const volScalarField& f
) const
{
    return averageOnTissueRegions(f, mesh_.V().field());
}


void Foam::tissueGraphRegions::writeFields() const
{
    nearestEdgeIndices_.write();
    nearestEdgeScalarIndices_.write();
}


// ************************************************************************* //
