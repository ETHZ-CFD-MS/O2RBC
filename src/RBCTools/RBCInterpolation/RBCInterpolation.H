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
    Foam::RBCInterpolation

Description
    Provides an interface for perform interpolation of fields between a
    collection of RBC objets and an Eulerian fixed grid.

SourceFiles
    RBCInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef RBCInterpolation_H
#define RBCInterpolation_H


#include "fvCFD.H"
#include "regIOobject.H"

#include "meshToMeshMoving.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class RBCCollection;

/*---------------------------------------------------------------------------*\
                        Class RBCInterpolation Declaration
\*---------------------------------------------------------------------------*/


class RBCInterpolation
:
    public regIOobject
{
private:
    // Private data

    //- Eulerian mesh
    const fvMesh& mesh_;

    //- List of RBCs
    RBCCollection& RBCCollection_;

    //- List of objects for interpolation from Eulerian to RBC mesh
    meshToMeshMoving eulerToRBCInterp_;

public:

    //- Runtime type information
    ClassName("RBCInterpolation");

    // Constructors

        // Construct from attributes
        RBCInterpolation
        (
            const fvMesh&,
            RBCCollection&
        );

    // Destructor
       ~RBCInterpolation();
    
    // Member functions

        // Access

            //- Constant access to the interpolation objects
            const meshToMeshMoving& eulerToRBCInterpolator() const
            {
                return eulerToRBCInterp_;
            }


        //- Update the interpolation objects
        void update();

        //- Interpolate a volVectorField defined on the Eulerian mesh to fields
        //  defined on RBCs
        void 
        interpolateEulerToRBC
        (
            const volScalarField& f_euler,
            const word& fieldName_RBC
        );

        //- Interpolate volVectorFields defined on RBCs to a field defined on
        // the Eulerian mesh.
        void 
        interpolateRBCToEuler
        (
            const word& fieldName_RBC,
            volScalarField& f_euler
        );
    
    //- Write
        bool writeData(Ostream&) const;

private:
        // Disable copy constructor
        RBCInterpolation
        (
            const RBCInterpolation&
        );

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

