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

Class
    Foam::deformableBodyGeometricState

Description
    Base class to describe the position and shape of a deformable body.

SourceFiles
    deformableBodyGeometricState.C

\*---------------------------------------------------------------------------*/

#ifndef deformableBodyGeometricState_H
#define deformableBodyGeometricState_H

#include "point.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "pointField.H"
#include "boundBox.H"

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class deformableBodyGeometricState;
Istream& operator>>(Istream&, deformableBodyGeometricState&);
Ostream& operator<<(Ostream&, const deformableBodyGeometricState&);


/*---------------------------------------------------------------------------*\
                Class deformableBodyGeometricState Declaration
\*---------------------------------------------------------------------------*/

class deformableBodyGeometricState
{
    // Private data

        //- Body center
        point center_;


public:

    //- Runtime type information
    TypeName("deformableBody");
        
    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            deformableBodyGeometricState,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            deformableBodyGeometricState,
            mesh,
            (
                const dictionary& dict,
                const polyMesh& mesh
            ),
            (dict, mesh)
        );


    // Constructors

        //- Construct from components
        deformableBodyGeometricState(const point& center);

        //- Construct from dictionary
        deformableBodyGeometricState(const dictionary&);

        //- Construct from Istream
        deformableBodyGeometricState(Istream&);

        //- Construct as copy
        deformableBodyGeometricState(const deformableBodyGeometricState&);

        //- Construct and return a clone
        virtual autoPtr<deformableBodyGeometricState> clone() const = 0;


    // Selectors

        //- Select constructed from dictionary
        static autoPtr<deformableBodyGeometricState> New(const dictionary&);

        //- Select constructed from dictionary and polyMesh
        static autoPtr<deformableBodyGeometricState> New
        (
            const dictionary&,
            const polyMesh&    
        );

    //- Destructor
    virtual ~deformableBodyGeometricState()
    {}


    // Member Functions

        // Access
            const point& center() const
            {
                return center_;
            }

            point& center()
            {
                return center_;
            }

            //- Return an approximate bounding box for the deformable body
            virtual boundBox bounds() const = 0;

        // Point transformation

            //- Transform points from the current geometric state to that given
            //  as argument
            virtual pointField transformPoints
            (
                const pointField& undisplacedPoints,
                const deformableBodyGeometricState& newState,
                const bool exactDiameter
            ) const = 0;
    

        // I/O
            virtual Istream& read(Istream&);
            virtual void write(Ostream&) const;


    // Member Operators

        void operator=(const deformableBodyGeometricState&);


    // IOstream Operators

        friend Istream& operator>>(Istream&, deformableBodyGeometricState&);
        friend Ostream& operator<<(Ostream&, const deformableBodyGeometricState&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
