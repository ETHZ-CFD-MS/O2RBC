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

#include "polygonalTube.H"

#include "geometricOps.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
Foam::polygonalTube::checkAndAdaptDiameters
(
    const dictionary& tubeOptions
)
{
    const bool useEffectiveDiameter = 
        readBool(tubeOptions.lookup("useEffectiveDiameter"));

    const scalar minimalDiameter = 
        tubeOptions.lookupOrDefault<scalar>
        (
            "minimalDiameter", 
            0.0
        );

    const scalar diameterFactor = 
        tubeOptions.lookupOrDefault<scalar>
        (
            "diameterFactor",
            1.0
        );

    if (innerDiameters_.size() != edge_.nSegments())
    {
        FatalErrorIn("polygonalTube::polygonalTube"
            "const polygonalEdge&, const scalarList&, "
            "const ellipseAxes&, const ellipseAxes&")
            << "The provided number of diameters (" << innerDiameters_.size()
            << ") does not match the number of edge segments ("
            << edge_.nSegments() << ")." << nl
            << abort(FatalError);
    }
    
    forAll(innerDiameters_, i)
    {
        innerDiameters_[i] = max(diameterFactor*innerDiameters_[i], minimalDiameter);
    }

    if (useEffectiveDiameter)
    {
        convertToEffectiveDiameter();
    }
}


void
Foam::polygonalTube::computeEllipseAxesList()
{
    axesList_.setSize(edge_.nPathVertices());

    axesList_.first() = startAxes();
    axesList_.last()  = endAxes();

    for (int i=1; i < edge_.nPathVertices() - 1; i++)
    {
        // find tangent vectors at path vertex
        vector v0 = edge_.tangentVectorFromPathVertex(i, false);
        vector v1 = edge_.tangentVectorFromPathVertex(i, true);

        // compute bisection vector and angle
        vector bisectV = Foam::bisectionVector(v0, v1);
        scalar bisectAngle = Foam::bisectionAngle(v0, v1);

        if (bisectAngle < constant::mathematical::pi/4.)
        {
            WarningIn("polygonalTube:computeEllipseAxesList()")
                << "The bisection angle at path vertex " << i 
                << " is too small." << nl
                << "Edge end points: " << edge_.startPoint() << ", "
                << edge_.endPoint() << nl
                << "Number of segments: " << edge_.nSegments() << endl;
        }

        // compute normal vector.
        // For stability reasons, use the bisection
        // vector to compute the normal vector. This yields better results when
        // v1 and v2 are nearly orthogonal.
        vector normalV = Foam::normalVector(v0, bisectV);

        scalar minorRadius = 0.25*(segmentOuterDiameters()[i-1] 
                                 + segmentOuterDiameters()[i]);
        scalar majorRadius = minorRadius/sin(bisectAngle);

        // construct ellipse axes
        axesList_[i] = ellipseAxes
                       (
                           normalV,
                           bisectV,
                           minorRadius,
                           majorRadius
                       );
    }
}


Foam::scalar
Foam::polygonalTube::effectiveDiameter() const
{
    const scalarList lengths = edge_.segmentLengths();

    // compute the scaled resistance of the sequence of circular segment
    scalar averageResistance = 0.0;

    forAll(innerDiameters_, i)
    {
        averageResistance += lengths[i]/pow(innerDiameters_[i],4);
    }
    averageResistance /= edge_.length();

    return pow(averageResistance, -0.25);
}


void
Foam::polygonalTube::convertToEffectiveDiameter()
{
    const scalar newDiameter = effectiveDiameter();
    forAll(innerDiameters_, i)
    {
        innerDiameters_[i] = newDiameter;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polygonalTube::polygonalTube
(
    const polygonalEdge& edge,
    const scalarList& innerDiameters,
    const ellipseAxes& startAxes,
    const ellipseAxes& endAxes,
    const dictionary& tubeOptions
)
:
    circularTube(startAxes, endAxes),
    edge_(edge),
    innerDiameters_(innerDiameters),
    outerInnerDiameterRatio_
    (
        readScalar(tubeOptions.lookup("outerInnerDiameterRatio"))
    ),
    axesList_()
{
    checkAndAdaptDiameters(tubeOptions);
    computeEllipseAxesList();
}


Foam::polygonalTube::polygonalTube
(
    const polygonalEdge& edge,
    const scalarList& innerDiameters,
    const dictionary& tubeOptions
)
:
    circularTube
    (
        ellipseAxes
        (
            Foam::normalVector(edge.tangentVector(0)),
            Foam::normalVector(edge.tangentVector(0), 
                               Foam::normalVector(edge.tangentVector(0))),
            0.5*innerDiameters[0],
            0.5*innerDiameters[0]
        ),
        ellipseAxes
        (
            Foam::normalVector(edge.tangentVector(edge.nSegments()-1)),
            Foam::normalVector(edge.tangentVector(edge.nSegments()-1), 
                               Foam::normalVector(edge.tangentVector(edge.nSegments()-1))),
            0.5*innerDiameters[edge.nSegments()-1],
            0.5*innerDiameters[edge.nSegments()-1]
        )
    ),
    edge_(edge),
    innerDiameters_(innerDiameters),
    outerInnerDiameterRatio_
    (
        readScalar(tubeOptions.lookup("outerInnerDiameterRatio"))
    ),
    axesList_()
{
    checkAndAdaptDiameters(tubeOptions);
    computeEllipseAxesList();
}


Foam::polygonalTube::polygonalTube
(
    const polygonalEdge& edge,
    const scalar innerDiameter,
    const ellipseAxes& startAxes,
    const ellipseAxes& endAxes,
    const scalar outerInnerDiameterRatio
)
:
    circularTube(startAxes, endAxes),
    edge_(edge),
    innerDiameters_(edge_.nSegments(), innerDiameter),
    outerInnerDiameterRatio_(outerInnerDiameterRatio),
    axesList_()
{
    computeEllipseAxesList();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polygonalTube::~polygonalTube()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::scalarList
Foam::polygonalTube::segmentOuterDiameters() const
{
    scalarList outerDiameters(innerDiameters_.size());
    forAll(innerDiameters_, i)
    {
        outerDiameters[i] = outerInnerDiameterRatio_*innerDiameters_[i];
    }
    return outerDiameters;
}


const Foam::scalarList
Foam::polygonalTube::segmentInnerDiameters() const
{
    return innerDiameters_;
}


void 
Foam::polygonalTube::setStartAxes(const ellipseAxes& startAxes)
{
    circularTube::setStartAxes(startAxes);
    computeEllipseAxesList();
}


void 
Foam::polygonalTube::setEndAxes(const ellipseAxes& endAxes)
{
    circularTube::setEndAxes(endAxes);
    computeEllipseAxesList();
}


Foam::scalar
Foam::polygonalTube::meanOuterDiameter() const
{
    return meanInnerDiameter()*outerInnerDiameterRatio_;
}


Foam::scalar
Foam::polygonalTube::meanInnerDiameter() const
{
    // do a length-weighted averaged of the diameter
    const scalarList lengths = edge_.segmentLengths();

    scalar meanDiameter = 0.0;
    forAll(innerDiameters_, i)
    {
        meanDiameter += innerDiameters_[i]*lengths[i];
    }
    return meanDiameter/edge_.length();
}


Foam::scalar
Foam::polygonalTube::outerDiameter(const scalar s) const
{
    label segmentI = edge_.segmentIndex(s);
    return segmentOuterDiameters()[segmentI];
}


Foam::scalar
Foam::polygonalTube::innerDiameter(const scalar s) const
{
    label segmentI = edge_.segmentIndex(s);
    return innerDiameters_[segmentI];
}


Foam::ellipseAxes
Foam::polygonalTube::ellipseAxesBefore(const scalar s) const
{
    const label segmentI = edge_.segmentIndex(s);
    return axesList_[segmentI];
}


Foam::ellipseAxes
Foam::polygonalTube::ellipseAxesAfter(const scalar s) const
{
    const label segmentI = edge_.segmentIndex(s);
    return axesList_[segmentI + 1];
}


Foam::ellipseAxes 
Foam::polygonalTube::ellipseAxesAtSegmentStart(const label segmentI) const
{
    return axesList_[segmentI];
}


Foam::ellipseAxes 
Foam::polygonalTube::ellipseAxesAtSegmentEnd(const label segmentI) const
{
    return axesList_[segmentI + 1];
}


Foam::pointField 
Foam::polygonalTube::skeletonPoints() const
{
    // Number of edge points for edge skeleton
    const label nEdgePoints = 4*edge_.nSegments();

    pointField pf;

    scalar eLength = edge().length();

    // add the point which are on the graph edge
    for(int i=0 ; i <= nEdgePoints ; i++)
    {
        scalar s = eLength*i/nEdgePoints;
        pf.append(edge().pointPosition(s));
    }

    // add the ellipse axes end points for each path vertex
    forAll(axesList_, i)
    {
        const ellipseAxes& cAxes = axesList_[i];
        const point& cPoint = edge().points()[i];
        pf.append
        (
            cPoint + cAxes.minorAxis()*cAxes.minorRadius()
        );
        pf.append
        (
            cPoint - cAxes.minorAxis()*cAxes.minorRadius()
        );
        pf.append
        (
            cPoint + cAxes.majorAxis()*cAxes.majorRadius()
        );
        pf.append
        (
            cPoint - cAxes.majorAxis()*cAxes.majorRadius()
        );
    }

    return pf;
}


void
Foam::polygonalTube::accept(tubeVisitor& visitor) const
{
    visitor.visitPolygonalTube(*this);
}


void
Foam::polygonalTube::write(Ostream& os) const
{
    circularTube::write(os);
    os  << "Edge: " << edge_
        << "Inner diameters: " << innerDiameters_
        << "Diameter ratio: " << outerInnerDiameterRatio_ << endl;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream&
Foam::operator<<(Ostream& os, const polygonalTube& obj)
{
    obj.write(os);

    return os;
}


// ************************************************************************* //
