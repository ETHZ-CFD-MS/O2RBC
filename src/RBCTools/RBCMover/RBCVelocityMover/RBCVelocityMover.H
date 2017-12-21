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
    Foam::RBCVelocityMover

Description
    Provide an interface to move RBCs for simulations in an 
    Eulerian frame of reference, given the RBC velocity, linear density and 
    spacings. This version works for axisymmetric simulations with a single 
    capillary.
    

SourceFiles
    RBCVelocityMover.C

\*---------------------------------------------------------------------------*/

#ifndef RBCVelocityMover_H
#define RBCVelocityMover_H

#include "fvCFD.H"

#include "RBCMover.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class randomVariable;

/*---------------------------------------------------------------------------*\
                        Class RBCVelocityMover Declaration
\*---------------------------------------------------------------------------*/


class RBCVelocityMover
:
    public RBCMover
{
protected:
    // Protected data

        //- Setter for inlet values
        autoPtr<randomVariable> inletValuePtr_;

        //- Name of the field defined at the inlet
        word inletFieldName_;

        //- Linear density
        scalar linearDensity_;

        //- RBC velocity
        scalar RBCVelocity_;

        //- Field name of the prescribed initial value
        word initialFieldName_;

        //- Initial value for RBCs in the domain
        scalar initialValue_;

        scalar RBCRadius_;
        scalar RBCVolume_;
        scalar plasmaRadius_;

        //- Name of the sample mesh used to generate RBC meshes
        word sampleMeshName_;


    // Protected member functions

        //- Compute the inlet PO2 value
        scalar computeInletValue() const;


public:

    //- Runtime type information
    TypeName("RBCVelocityMover");

    
    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            RBCVelocityMover,
            meshRBCCollection,
            (
                const fvMesh& mesh,
                RBCCollection& RBCCollection
            ),
            (mesh, RBCCollection)
        );


    // Selectors

        //- Return a reference to the selected interpolation scheme
        static autoPtr<RBCVelocityMover> New
        (
            const fvMesh& mesh,
            RBCCollection& RBCCollection
        );
   

    // Constructors

        //- Construct from attributes
        RBCVelocityMover
        (
            const fvMesh& mesh,
            RBCCollection& RBCCollection
        );

    // Destructor
       ~RBCVelocityMover();
    
    // Member functions

        //- Set initial positions of RBC meshes
        void setInitialPositions();

        //- Move all RBCs using the velocity field U
        void moveAll();

        //- Compute the current RBC diameter
        virtual scalar RBCDiameter() const;

        //- Compute the current RBC radius
        virtual scalar RBCRadius() const;

        //- Setter for the linear density
            void setLinearDensity(const scalar linearDensity) 
            {
                linearDensity_ = linearDensity;
            }

        //- Setter for the RBC velocity
            void setRBCVelocity(const scalar RBCVelocity)
            {
                RBCVelocity_ = RBCVelocity;
            }

        //- Getter for the linear density
            scalar getLinearDensity() const 
            {
                return linearDensity_;
            }

        //- Getter for the RBC velocity
            scalar getRBCVelocity() const 
            {
                return RBCVelocity_;
            }

        //- Write
            bool writeData(Ostream&) const;

private:
        //- Disable copy constructor
        RBCVelocityMover
        (
            const RBCVelocityMover&
        );

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

