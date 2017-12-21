/*---------------------------------------------------------------------------*\

Class
    Foam::triSurfaceDistanceFunction

Description
    A class to evaluate volScalarFields based on the distance to a triSurface.

SourceFiles
    triSurfaceDistanceFunction.C

Author
    Adrien Luecker, Jul 2013


\*---------------------------------------------------------------------------*/

#ifndef triSurfaceDistanceFunction_H
#define triSurfaceDistanceFunction_H

#include "polyMesh.H"
#include "cellSet.H"
#include "topoSet.H"
#include "topoSetSource.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class triSurfaceSearch;
class triSurface;

/*---------------------------------------------------------------------------*\
               Class triSurfaceDistanceFunction Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceDistanceFunction
{
    // Private data

        //- Mesh
        const polyMesh& mesh_; 
    
        //- Name of surface file
        const fileName surfName_;

        //- includes into the computation cells with distance from cellCentre
        //  to surface less than threshold
        const scalar threshold_;

        //- Function type
        //  Default type: Gaussian
        const word functionType_;

        //- scalar parameters needed to evaluate the function
        const List<scalar> params_;

        //- points outside the surface
        const pointField outsidePoints_;

        //- Cells which are near the surface
        cellSet nearSurfaceCells_;

        //- Cells which are inside the surface
        cellSet insideSurfaceCells_;

        //- Cells which are outside the surface
        cellSet outsideSurfaceCells_;

        //- Distance from cell centers to the surface
        scalarField distanceToSurface_;
        
        //- triSurface to use for distance calculation
        const triSurface* surfPtr_;
        
        //- search engine on surface.
        const triSurfaceSearch* querySurfPtr_;

        //- whether I allocated above surface ptrs or whether they are
        //  external.
        bool IOwnPtrs_;

    // Private Member Functions

        //- Build cellSet nearSurfaceCells_
        void buildNearSurfaceCellSet();

        //- Compute distance to surface
        void computeDistanceToSurface();
        
public:
    // Constructors

        //- Construct from components
        triSurfaceDistanceFunction
        (
            const polyMesh& mesh,
            const fileName& surfName,
            const scalar threshold,
            const word functionType,
            const List<scalar> params,
            const pointField& outsidePoints
        );
        
        //- Construct from dictionary
        triSurfaceDistanceFunction
        (
            const polyMesh& mesh,
            const IOdictionary& dict
        );

    //- Destructor
    virtual ~triSurfaceDistanceFunction();
    


    // Member Functions

        virtual void evaluate
        (
            volScalarField& f
        );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
