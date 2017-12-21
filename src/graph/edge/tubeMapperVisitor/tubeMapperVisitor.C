/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "tubeMapperVisitor.H"

#include "geometricOps.H"
#include "dimensionSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tubeMapperVisitor, 0);
    const word Foam::tubeMapperVisitor::tubeMeshPrefix_ = "tube";
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word
Foam::tubeMapperVisitor::tubeMeshName(const label edgeIndex)
{
    std::stringstream sstm;
    sstm << tubeMeshPrefix_ << edgeIndex;
    return sstm.str();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
Foam::tubeMapperVisitor::computeCylinderData()
{
    boundBox bb = cylinderMesh_.bounds();

    cylinderAxisPoints_.first()  = point(bb.min().x(), bb.midpoint().y(), bb.midpoint().z());
    cylinderAxisPoints_.second() = point(bb.max().x(), bb.midpoint().y(), bb.midpoint().z());
    cylinderRadius_ = 0.5*bb.span().y();

    if (debug)
    {
        Info<< "tubeMapperVisitor::computeCylinderData(): " << nl
            << "Left axis point: " << cylinderAxisPoints_.first() << nl
            << "Right axis point: " << cylinderAxisPoints_.second() << nl
            << "Cylinder radius: " << cylinderRadius_ << endl;
    }
}


void
Foam::tubeMapperVisitor::computeGeometricTransformations
(
    const polygonalTube& tube
)
{
    rotationToAxisList_.setSize(tube.edge().nSegments());

    forAll(rotationToAxisList_, segmentI)
    {
        rotationToAxisList_[segmentI] = 
            rotationTensorFromSegmentToAxis(tube, segmentI);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tubeMapperVisitor::tubeMapperVisitor
(
    const IOobject& io
)
:
    tubeVisitor(),
    cylinderMesh_
    (
        IOobject
        (
            io.name(),
            io.instance(),
            io.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    ),
    transformedPoints_(cylinderMesh_.points()),
    cylinderAxisPoints_(),
    cylinderRadius_(),
    rotationToAxisList_()
{
    computeCylinderData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tubeMapperVisitor::~tubeMapperVisitor()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyMesh> 
Foam::tubeMapperVisitor::mappedMesh(const label tubeIndex) const
{
    return autoPtr<polyMesh>
    (
        new polyMesh
        (
            IOobject
            (
                tubeMeshName(tubeIndex),
                cylinderMesh_.time().timeName(),
                cylinderMesh_.time(),
                IOobject::NO_READ
            ),
            xferCopy(transformedPoints_),
            xferCopy(cylinderMesh_.faces()),
            xferCopy(cylinderMesh_.faceOwner()),
            xferCopy(cylinderMesh_.faceNeighbour()),
            false // no parallel comms
        )
    );
}


Foam::autoPtr<volScalarField>
Foam::tubeMapperVisitor::normalizedCylinderAxisCoord() const
{
    const dimensionedScalar xLeft( "xLeft", dimLength, cylinderAxisPoints_.first().x());
    const dimensionedScalar xRight("xRight", dimLength, cylinderAxisPoints_.second().x());
    
    autoPtr<volScalarField> pCoordField
    (
        new volScalarField
        (
            IOobject
            (
                "xAxisCoord",
                cylinderMesh_.time().timeName(),
                cylinderMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            cylinderMesh_,
            dimless
        )
    );

    pCoordField->primitiveFieldRef() = 
        (cylinderMesh_.C().component(0) - xLeft)/(xRight - xLeft);
    forAll(cylinderMesh_.boundary(), patchI)
    {
        const polyPatch& patch = cylinderMesh_.boundaryMesh()[patchI];
        const vectorField::subField faceCentres = patch.faceCentres();

        forAll(patch, i)
        {
            pCoordField->boundaryFieldRef()[patchI][i] = 
                (faceCentres[i][0] - xLeft.value())/(xRight - xLeft).value();

        }
    }

    return pCoordField;
}


void 
Foam::tubeMapperVisitor::visitPolygonalTube(const polygonalTube& tube)
{
    computeGeometricTransformations(tube);

    for (int segmentI = 0; segmentI < tube.edge().nSegments(); segmentI++)
    {
        if (debug)
        {
            Info<< "tubeMapperVisitor::visit(const polygonalTube&): "
                << " dealing with segment " << segmentI << endl;
        }

        // get indices of cylinder points in current segment
        labelList pointIDs = cylinderPointIDsOnSegment(tube, segmentI);

        // transform all the points in the given segment
        forAll(pointIDs, i)
        {
            label pI = pointIDs[i];

            transformedPoints_[pI] = 
                transformPoint
                (
                    cylinderMesh_.points()[pI],
                    tube,
                    segmentI
                );
        }
    }
}


Foam::point
Foam::tubeMapperVisitor::transformPoint
(
    const point& p,
    const polygonalTube& tube,
    const label segmentI
) const
{
    const ellipseAxes axesStart = tube.ellipseAxesAtSegmentStart(segmentI);
    const ellipseAxes axesEnd   = tube.ellipseAxesAtSegmentEnd  (segmentI);

    // build transformation tensor from segment axis to x-axis
    // const tensor rotationToAxis = rotationTensorFromSegmentToAxis(tube, segmentI);
    const tensor rotationToAxis = rotationToAxisList_[segmentI];

    // build shear tensors that correspond to respective ellipseAxes
    const tensor shearLeft = shearTensorOnCylinderAxis
    (
        axesStart,
        rotationToAxis
    );
    const tensor shearRight = shearTensorOnCylinderAxis
    (
        axesEnd,
        rotationToAxis
    );

    const point cylinderAxisLeft  = cylinderAxisPoints_.first();
    const point cylinderAxisRight = cylinderAxisPoints_.second();
    const scalar cylinderLength = cylinderAxisRight.x() - cylinderAxisLeft.x();

    const scalar edgeLength    = tube.edge().length();
    const scalar segmentLength = tube.edge().segmentLengths()[segmentI];

    // construct scaled projection of p on the left/right sides of the
    // scaled cylinder
    point pLeft  = point
                   (  
                       0.0,
                       p.y()*axesStart.minorRadius()/cylinderRadius_,
                       p.z()*axesStart.minorRadius()/cylinderRadius_
                   );
    point pRight = point
                   (  
                       segmentLength,
                       p.y()*axesEnd.minorRadius()/cylinderRadius_,
                       p.z()*axesEnd.minorRadius()/cylinderRadius_
                   );
    // shear projections so that they are on the left/right side of the elliptic
    // cylinder
    pLeft  = shearLeft  & pLeft;
    pRight = shearRight & pRight;

    // weighting factor for convex combination
    scalarList vertexCoords = tube.edge().pathVertexCoords();
    scalar lambdaInEdge    = (p.x() - cylinderAxisLeft.x())
                             *edgeLength/cylinderLength;
    scalar lambdaInSegment = (lambdaInEdge - vertexCoords[segmentI])
                             /segmentLength;

    // compute convex combination of pLeft and pRight
    point mappedPoint = (1.0 - lambdaInSegment)*pLeft + lambdaInSegment*pRight;

    // mapped the point to the tube mesh
    mappedPoint = tube.edge().points()[segmentI]
                 + det(rotationToAxis)*(inv(rotationToAxis) & mappedPoint);

    return mappedPoint;
}


Foam::tensor
Foam::tubeMapperVisitor::rotationTensorFromSegmentToAxis
(
    const polygonalTube& tube,
    const label segmentI
) const
{
    vector segmentAxis = tube.edge().tangentVector(segmentI);
    vector axisNormal  = tube.ellipseAxesAtSegmentStart(segmentI).minorAxis();

    if (segmentI > 0)
    {
        vector previousAxis = tube.edge().tangentVector(segmentI - 1);

        return rotationTensorFromSegmentToAxis(tube, segmentI - 1)
               & Foam::rotationTensor(segmentAxis, previousAxis);
    }
    else
    {
        return Foam::rotationTensor(Pair<vector>(segmentAxis, axisNormal),
                                    Pair<vector>(vector(1,0,0), vector(0,0,1))); 
    }
}


Foam::tensor
Foam::tubeMapperVisitor::shearTensorOnCylinderAxis
(
    const ellipseAxes& axes,
    const tensor& T
) const
{
    vector transformedMajorAxis = T & axes.majorAxis();
    vector projectedToYZPlane = transformedMajorAxis;
    projectedToYZPlane.x() = 0.0;

    return shearTensor(projectedToYZPlane, transformedMajorAxis);
}


Foam::labelList 
Foam::tubeMapperVisitor::cylinderPointIDsOnSegment
(
    const polygonalTube& tube,
    const label segmentI
) const
{
    labelList pointIDs;

    const scalar cylinderLength = cylinderAxisPoints_.second().x()
                          - cylinderAxisPoints_.first().x();
    const scalar edgeLength = tube.edge().length();

    forAll(cylinderMesh_.points(), pI)
    {
        const point& p = cylinderMesh_.points()[pI];
        const scalar pointS = (p.x() - cylinderAxisPoints_.first().x())
                              *edgeLength/cylinderLength;

        // optimize this call
        const label pointSegmentI = tube.edge().segmentIndex(pointS);

        if (pointSegmentI == segmentI)
        {
            pointIDs.append(pI);
        }
    }

    return pointIDs;
}

// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
