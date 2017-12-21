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
    Foam::RBCMover

Description
    Abstract class for moving RBCs through a fixed mesh.

SourceFiles
    RBCMover.C

\*---------------------------------------------------------------------------*/

#ifndef RBCMover_H
#define RBCMover_H

#include "fvCFD.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class RBCCollection;


/*---------------------------------------------------------------------------*\
                        Class RBCMover Declaration
\*---------------------------------------------------------------------------*/


class RBCMover
:
    public IOdictionary
{
protected:
    // Protected data

    //- Eulerian mesh
    const fvMesh& mesh_;

    //- List of RBCs
    RBCCollection& RBCCollection_;


    // Protected member functions
    
        //- Compute PO2 based on the field name and its value
        scalar computePO2FromFieldNameAndValue
        (
            const scalar, 
            const word& fieldName
        ) const;


public:

    //- Runtime type information
    TypeName("RBCMover");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            RBCMover,
            meshRBCCollection,
            (
                const fvMesh& mesh,
                RBCCollection& RBCCollection
            ),
            (mesh, RBCCollection)
        );


    // Selectors

        //- Return a reference to the selected interpolation scheme
        static autoPtr<RBCMover> New
        (
            const fvMesh& mesh,
            RBCCollection& RBCCollection
        );
   

    // Constructors

        //- Construct from attributes
        RBCMover
        (
            const fvMesh& mesh,
            RBCCollection& RBCCollection
        );

    // Destructor
       ~RBCMover();
    
    // Member functions

        //- Set initial positions of RBC meshes
        virtual void setInitialPositions() = 0;

        //- Move the RBCs
        virtual void moveAll() = 0;

        //- Write
            bool writeData(Ostream&) const;

private:
        //- Disable copy constructor
        RBCMover(const RBCMover&);

        //- Disallow default bitwise assignment
        void operator=(const RBCMover&);

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

