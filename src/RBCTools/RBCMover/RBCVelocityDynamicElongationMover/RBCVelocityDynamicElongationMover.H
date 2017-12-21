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
    Foam::RBCVelocityDynamicElongationMover

Description
    Provide an interface to move RBCs for simulations in an 
    Eulerian frame of reference, given the RBC velocity, linear density and 
    spacings. This version works for axisymmetric simulations with a single 
    capillary.
    

SourceFiles
    RBCVelocityDynamicElongationMover.C

\*---------------------------------------------------------------------------*/

#ifndef RBCVelocityDynamicElongationMover_H
#define RBCVelocityDynamicElongationMover_H

#include "RBCVelocityMover.H"

#include "interpolationTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class RBCVelocityDynamicElongationMover Declaration
\*---------------------------------------------------------------------------*/


class RBCVelocityDynamicElongationMover
:
    public RBCVelocityMover
{
private:
    // Private data
    
        // Tube hematocrit
        scalar tubeHematocrit_;

        //  Reference RBC velocity
        scalar RBCVelocityReference_;

        //- Relation between PO2 and the RBC radius
        interpolationTable<scalar> RBCRadiusTable_;

        //- Relation between PO2 and the RBC velocity factor
        interpolationTable<scalar> RBCVelocityFactorTable_;

        //- Name of the cell zone for averaging PO2
        word PO2AveragingZone_;


    // Private member functions

        //- Return the PO2 average used to evaluate the RBC radius and the RBC
        //  velocity factor
        scalar PO2Average() const;


public:

    //- Runtime type information
    ClassName("RBCVelocityDynamicElongationMover");

    //- Static data members
    static const word plasmaZoneName;

    
    // Constructors

        //- Construct from attributes
        RBCVelocityDynamicElongationMover
        (
            const fvMesh& mesh,
            RBCCollection& RBCCollection
        );

    // Destructor
       ~RBCVelocityDynamicElongationMover();
    
    // Member functions

        //- Move all RBCs using the velocity field U
        void moveAll();

        //- Compute the current RBC radius
        scalar RBCRadius() const;

private:
        //- Disable copy constructor
        RBCVelocityDynamicElongationMover
        (
            const RBCVelocityDynamicElongationMover&
        );
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

