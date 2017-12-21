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

#include "polygonalEdge.H"

#include "Pair.H"
#include "Field.H"
#include "ListOps.H"
#include "ops.H"

#include "polygonalTube.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
Foam::polygonalEdge::computeLength()
{
    length_ = Foam::sum(segmentLengths());
}


void
Foam::polygonalEdge::computePathVertexCoords()
{
    pathVertexCoords_.setSize(nPathVertices());

    scalarList lengths = segmentLengths();

    pathVertexCoords_[0] = 0.0;

    for (int i = 1; i < pathVertexCoords_.size(); i++)
    {
        pathVertexCoords_[i] = pathVertexCoords_[i-1] + lengths[i-1];
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polygonalEdge::polygonalEdge(const List<point>& points)
:
    geometricEdge(),
    points_(points)
{
    points_.checkSize(1);

    startPoint() = points_[0];
    endPoint() = points_.last();

    computeLength();
    computePathVertexCoords();
}


Foam::polygonalEdge::polygonalEdge
(
    const point& start,
    const point& end
)
:
    geometricEdge(start, end),
    points_(2)
{
    points_[0] = start;
    points_[1] = end;

    computeLength();
    computePathVertexCoords();
}


Foam::polygonalEdge::polygonalEdge(const Pair<point>& pp)
:
    geometricEdge(pp),
    points_(2)
{
    points_[0] = pp.first();
    points_[1] = pp.second();

    computeLength();
    computePathVertexCoords();
}


Foam::polygonalEdge::polygonalEdge(Istream& is)
:
    geometricEdge(),
    points_(is)
{
    points_.checkSize(1);

    startPoint() = points_[0];
    endPoint() = points_.last();

    computeLength();
    computePathVertexCoords();
}


Foam::polygonalEdge::polygonalEdge(const polygonalEdge& obj)
:
    geometricEdge(obj),
    points_(obj.points_)
{
    computePathVertexCoords();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polygonalEdge::~polygonalEdge()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::point
Foam::polygonalEdge::pointPosition(const scalar s) const
{
    const scalar clampedS = checkAndClampCurvilinearCoord(s);

    const label segmentI = segmentIndex(clampedS);
    const scalar lambda = (s - pathVertexCoords()[segmentI])/segmentLengths()[segmentI];

    const point p1 = points_[segmentI];
    const point p2 = points_[segmentI + 1];

    return (1.0 - lambda)*p1 + lambda*p2;
}


Foam::vector
Foam::polygonalEdge::tangentVector(const scalar s) const
{
    const Pair<point> pp = segmentEndPoints(s);
    const vector edgeV = pp.second() - pp.first();

    if (edgeV == vector::zero)
    {
        FatalErrorIn
        (
            "Foam::polygonalEdge::tangentVector(const scalar)"
        ) << "The intermediate points " << pp.first() 
          << " and " << pp.second() << " have the "
          << "same positions. The edge tangent vector is undefined." << nl
          << abort(FatalError);
    }

    return edgeV/mag(edgeV);
}


Foam::vector
Foam::polygonalEdge::tangentVector(const label segmentI) const
{
    return tangentVectorFromPathVertex(segmentI, true);
}



Foam::Pair<Foam::point>
Foam::polygonalEdge::segmentEndPoints(const scalar s) const
{
    const scalar clampedS = checkAndClampCurvilinearCoord(s);
    const label segmentI = segmentIndex(clampedS);

    return Pair<point>(points_[segmentI], points_[segmentI + 1]);
}


Foam::label
Foam::polygonalEdge::segmentIndex(const scalar s) const
{
    const scalar clampedS = checkAndClampCurvilinearCoord(s);
    const scalarList pathCoords = pathVertexCoords();

    if (clampedS < pathCoords[0])
    {
        return 0;
    }
    else if (clampedS >= pathCoords.last())
    {
        return (nSegments() - 1);
    }
    else
    {
        return findLower(pathCoords, s, 0, Foam::lessEqOp<scalar>());
    }

    return -1;
}


Foam::scalarList
Foam::polygonalEdge::segmentLengths() const
{
    scalarList segmentLengths(nSegments());

    forAll(segmentLengths, i)
    {
        point p1 = points_[i];
        point p2 = points_[i+1];
        segmentLengths[i] = mag(p2 - p1);
    }

    return segmentLengths;
}


const Foam::scalarList&
Foam::polygonalEdge::pathVertexCoords() const
{
    return pathVertexCoords_;
}


Foam::vector
Foam::polygonalEdge::tangentVectorFromPathVertex
(
    const label pathVertexI, 
    const bool positiveDirection
) const
{
    if (pathVertexI == 0 && !(positiveDirection))
    {
        FatalErrorIn("tangentVectorFromPathVertex(const label, const bool")
            << "The tangent vector from the first path vertex in the negative "
            << "direction is undefined." << nl
            << abort(FatalError);
    }
    else if (pathVertexI == nPathVertices() && positiveDirection)
    {
        FatalErrorIn("tangentVectorFromPathVertex(const label, const bool")
            << "The tangent vector from the last path vertex in the positive "
            << "direction is undefined." << nl
            << abort(FatalError);
    }

    vector tangentV;

    // find desired curvilinear coordinate
    scalar s;

    if (positiveDirection)
    {
        s = 0.5*(pathVertexCoords_[pathVertexI    ] 
               + pathVertexCoords_[pathVertexI + 1]);
        tangentV = tangentVector(s);
    }
    else
    {
        s = 0.5*(pathVertexCoords_[pathVertexI - 1] 
               + pathVertexCoords_[pathVertexI    ]);
        tangentV = -tangentVector(s);
    }

    return tangentV;
}


Foam::autoPtr<Foam::circularTube>
Foam::polygonalEdge::makeCircularTube
(
    const scalarList& segmentDiameters,
    const dictionary& tubeOptions
) const
{
    autoPtr<circularTube> pTube
    (
        new polygonalTube
        (
            *this,
            segmentDiameters,
            tubeOptions
        )
    );

    return pTube;
}


void Foam::polygonalEdge::write(Ostream& os) const
{
    geometricEdge::write(os);
    os << "Path vertices: " << points_ << endl;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::polygonalEdge::operator=(const polygonalEdge& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::polygonalEdge::operator=(const Foam::polygonalEdge&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    this->geometricEdge::operator=(rhs);
    points_ = rhs.points();
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream&
Foam::operator<<(Ostream& os, const polygonalEdge& obj)
{
    obj.write(os);

    return os;
}


// ************************************************************************* //
