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
    Foam::RBCPath

Description
    Provide an interface to move RBCs for simulations in an 
    Eulerian frame of reference. The motion is defined by externally defined
    RBC paths.
    

SourceFiles
    RBCPathMover.C

\*---------------------------------------------------------------------------*/

#ifndef RBCPathMover_H
#define RBCPathMover_H

#include "RBCMover.H"

#include "fvCFD.H"
#include "regIOobject.H"

#include "RBCPathCollection.H"
#include "graphCoordinateInterpolation.H"
#include "deformableBodyGeometricState.H"
#include "RBCProperty.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
// class RBC;
class RBCPath;

class circularTubeGraph;
class diameterFunction;
class graphInletValue;


/*---------------------------------------------------------------------------*\
                        Class RBCPathMover Declaration
\*---------------------------------------------------------------------------*/


class RBCPathMover
:
    public RBCMover
{
private:
    // Private data

        //- RBC paths
        RBCPathCollection RBCPaths_;

        //- Reference to graph
        const circularTubeGraph& graph_;

        //- Interpolation object for the RBC paths
        graphCoordinateInterpolation graphInterpolator_;

        //- Name of the sample meshes for the RBCs
        RBCProperty<word> sampleMeshNames_;

        //- Function used to compute the RBC diameter
        autoPtr<diameterFunction> RBCDiameterFunctionPtr_;

        //- Setter for inlet values
        autoPtr<graphInletValue> inletValuePtr_;

        //- Field name of the prescribed initial value
        word initialFieldName_;

        //- Initial value for RBCs in the domain
        scalar initialValue_;

        //- RBC volume
        scalar RBCVolume_;


    // Private member functions

        //- Compute the current position and direction on a RBC Path
        autoPtr<deformableBodyGeometricState> computeRBCGeometricState
        (
            const RBCPath&
        ) const;

        //- Compute the inlet PO2 value for a given RBC path
        scalar computeInletValue(const RBCPath&) const;

        //- Return the inlet value associated to the edge where the RBC path
        //  currently resides
        scalar currentInletValue(const RBCPath&) const;

        //- Normalize the curvilinear coordinates of the RBC paths from 0 to 1.
        void normalizePathCoordinates();

        //- Convert the path edge indices from external to internal
        void convertPathEdgeIndices();

        //- Disallow default bitwise copy construct
        RBCPathMover(const RBCPathMover&);

        //- Disallow default bitwise assignment
        void operator=(const RBCPathMover&);


public:

    //- Runtime type information
    ClassName("RBCPathMover");

    // Constructors

        //- Construct from attributes and reads graph from database
        RBCPathMover
        (
            const fvMesh& mesh,
            RBCCollection& RBCCollection
        );

    // Destructor
       ~RBCPathMover();
    
    // Member functions

        // Access

            const RBCPathCollection& RBCPaths() const
            {
                return RBCPaths_;
            }


        //- Set initial positions of RBC meshes
        void setInitialPositions();

        //- Move all RBCs
        void moveAll();

        //- Graph coordinate of RBC
        graphCoordinateDirected graphCoordinate(const label pathI) const;

        //- Write
        Ostream& write(Ostream&) const;


    // IOstream operator
        
        friend Ostream& operator<<(Ostream&, const RBCPathMover&);

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


