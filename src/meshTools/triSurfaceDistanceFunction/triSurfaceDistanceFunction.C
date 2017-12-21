#include "triSurfaceDistanceFunction.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "cellSet.H"
#include "surfaceToCell.H"
#include "topoSetSource.H"
#include "IOobject.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::triSurfaceDistanceFunction::buildNearSurfaceCellSet()
{
    surfaceToCell nearSurfaceSetSource
    (
       mesh_,
       surfName_,
       *surfPtr_,
       *querySurfPtr_,
       outsidePoints_,
       false, // include cut
       false, // include inside
       false, // include outside
       false, // include surface orientation
       threshold_,
       -100
    );

    nearSurfaceSetSource.applyToSet(topoSetSource::ADD, nearSurfaceCells_);

    surfaceToCell insideSurfaceSetSource
    (
       mesh_,
       surfName_,
       *surfPtr_,
       *querySurfPtr_,
       outsidePoints_,
       true,
       true, // include inside
       false,
       false, // include surface orientation
       -1,
       -100
    );

    insideSurfaceSetSource.applyToSet(topoSetSource::ADD, insideSurfaceCells_);

    surfaceToCell outsideSurfaceSetSource
    (
       mesh_,
       surfName_,
       *surfPtr_,
       *querySurfPtr_,
       outsidePoints_,
       true,
       false,
       true, // include outside
       false, // include surface orientation
       -1,
       -100
    );

    outsideSurfaceSetSource.applyToSet(topoSetSource::ADD, outsideSurfaceCells_);
}

void Foam::triSurfaceDistanceFunction::computeDistanceToSurface()
{
    // initialize all values to a negative number (only cells within a
    // distance less than threshold will have positive values)
    distanceToSurface_ = 100;

    const pointField& ctrs = mesh_.cellCentres();
    
    // Box dimensions to search in octree.
    const vector span(threshold_, threshold_, threshold_);

    forAll(ctrs, cellI)
    {
        const point& c = ctrs[cellI];
        // WARNING: since STL surfaces are three dimensional, in case of
        // extruded surfaces (for 2D simulation), it is possible that the
        // hit point is on the plane which is orthogonal to the extrusion
        // direction.
        // Accordingly, the extrusion length of the surface has to be chosen
        // large enough (greater than 2*threshold_ to be on the safe side).
        pointIndexHit inter = querySurfPtr_->nearest(c, span);

        if (inter.hit())
        {
            // compute the distance using just x- and y-dimensions
            const scalar dist = sqrt(sqr(inter.hitPoint().x() - c.x())
                                   + sqr(inter.hitPoint().y() - c.y()));
            distanceToSurface_[cellI] = dist;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceDistanceFunction::triSurfaceDistanceFunction
(
    const polyMesh& mesh,
    const fileName& surfName,
    const scalar threshold,
    const word functionType,
    const List<scalar> params,
    const pointField& outsidePoints
)
:
    mesh_(mesh),
    surfName_(surfName),
    threshold_(threshold),
    functionType_(functionType),
    params_(params),
    outsidePoints_(outsidePoints),
    nearSurfaceCells_(mesh_, "nearSurfaceCells", IOobject::NO_READ, IOobject::AUTO_WRITE),
    insideSurfaceCells_(mesh_, "insideSurfaceCells", IOobject::NO_READ, IOobject::AUTO_WRITE),
    outsideSurfaceCells_(mesh_, "outsideSurfaceCells", IOobject::NO_READ, IOobject::AUTO_WRITE),
    distanceToSurface_(mesh_.cellCentres().size()),
    surfPtr_(new triSurface(surfName_)),
    querySurfPtr_(new triSurfaceSearch(*surfPtr_)),
    IOwnPtrs_(true)
{
    buildNearSurfaceCellSet();
    computeDistanceToSurface();
}

Foam::triSurfaceDistanceFunction::triSurfaceDistanceFunction
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    mesh_(mesh),
    surfName_(dict.lookup("file")),
    threshold_(readScalar(dict.lookup("threshold"))),
    functionType_(dict.lookup("functionType")),
    params_(dict.lookup("params")),
    outsidePoints_(dict.lookup("outsidePoints")),
    nearSurfaceCells_(mesh_, "nearSurfaceCells", IOobject::NO_READ, IOobject::AUTO_WRITE),
    insideSurfaceCells_(mesh_, "insideSurfaceCells", IOobject::NO_READ, IOobject::AUTO_WRITE),
    outsideSurfaceCells_(mesh_, "outsideSurfaceCells", IOobject::NO_READ, IOobject::AUTO_WRITE),
    distanceToSurface_(mesh_.cellCentres().size()),
    surfPtr_(new triSurface(surfName_)),
    querySurfPtr_(new triSurfaceSearch(*surfPtr_)),
    IOwnPtrs_(true)
{
    buildNearSurfaceCellSet();
    computeDistanceToSurface();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceDistanceFunction::~triSurfaceDistanceFunction()
{
    if (IOwnPtrs_)
    {
        if (surfPtr_)
        {
            delete surfPtr_;
        }
        if (querySurfPtr_)
        {
            delete querySurfPtr_;
        }
    }
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurfaceDistanceFunction::evaluate
(
    volScalarField& f
)
{
    Field<scalar>& f_internal = f.primitiveFieldRef();

    if (functionType_ == "Gaussian")
    {
        scalar baseValue = params_[0]; // value away from surface
        scalar delta     = params_[1]; // height of the Gaussian "perturbation"
        scalar sigma     = params_[2]; // std. dev. of Gaussian curve

        f_internal = baseValue + delta*exp(-sqr(distanceToSurface_/sigma));
    }
    else if (functionType_ == "Gaussian_exterior")
    {
        scalar delta     = params_[0]; // height of the Gaussian "perturbation"
        scalar sigma     = params_[1]; // std. dev. of Gaussian curve

        const pointField& ctrs = mesh_.cellCentres();

        forAllConstIter(topoSet, outsideSurfaceCells_, iter)
        {
            const label cellI = iter.key();

            if (distanceToSurface_[cellI] <= threshold_
                && ctrs[cellI].y() <= 1.8e-6)
            {
                f_internal[cellI] += delta*exp(-sqr(distanceToSurface_[cellI]/sigma));
            }
        }
    }
}






