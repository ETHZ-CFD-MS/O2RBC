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
    Foam::axisymmetricBodyGeometricState

Description
    Describes the position of an axisymmetric body with given center, axis,
    diameter and length.

    The axis vector is not assumed to be normalized.

SourceFiles
    axisymmetricBodyGeometricState.C

\*---------------------------------------------------------------------------*/

#ifndef axisymmetricBodyGeometricState_H
#define axisymmetricBodyGeometricState_H

#include "deformableBodyGeometricState.H"

#include "autoPtr.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class axisymmetricBodyGeometricState;
Istream& operator>>(Istream&, axisymmetricBodyGeometricState&);
Ostream& operator<<(Ostream&, const axisymmetricBodyGeometricState&);


/*---------------------------------------------------------------------------*\
                Class axisymmetricBodyGeometricState Declaration
\*---------------------------------------------------------------------------*/

class axisymmetricBodyGeometricState
:
    public deformableBodyGeometricState
{
    // Private data

        //- Body axis
        vector axis_;

        //- Body diameter
        scalar diameter_;

        //- Body length on the centerline
        scalar length_;

public:

    //- Runtime type information
    TypeName("axisymmetricBody");


    // Constructors

        //- Construct from components
        axisymmetricBodyGeometricState
        (
            const point& center,
            const vector& axis,
            const scalar diameter,
            const scalar length
        );

        //- Construct from dictionary
        axisymmetricBodyGeometricState(const dictionary&);

        //- Construct from dictionary and polyMesh
        axisymmetricBodyGeometricState
        (
            const dictionary&,
            const polyMesh&
        );

        //- Construct from Istream
        axisymmetricBodyGeometricState(Istream&);

        //- Construct as copy
        axisymmetricBodyGeometricState(const axisymmetricBodyGeometricState&);

        //- Construct and return a clone
        virtual autoPtr<deformableBodyGeometricState> clone() const
        {
            return autoPtr<deformableBodyGeometricState>
            (
                new axisymmetricBodyGeometricState(*this)
            );
        }


    //- Destructor
    virtual ~axisymmetricBodyGeometricState()
    {}


    // Member Functions

        // Access

            const vector& axis() const
            {
                return axis_;
            }

            vector& axis()
            {
                return axis_;
            }

            const scalar& diameter() const
            {
                return diameter_;
            }

            scalar& diameter()
            {
                return diameter_;
            }

            const scalar& length() const
            {
                return length_;
            }

            scalar& length()
            {
                return length_;
            }

            //- Return an approximate bounding box of the axisymmetric body.
            //  This is done by computing the bounding box of a cylinder with
            //  the given axis, length, center and diameter. The center point is
            //  assumed to be on the middle of the axis, at distance length_/2
            //  of the axis end points.
            boundBox bounds() const;

        // Point transformation

            //- Transform points from the current geometric state to that given
            //  as argument. The transformation conserves the body volume. 
            //  If exactDiameter is true, the transformation is set so that the
            //  resulting diameter exactly corresponds to that of newState.
            //  If it is false, the resulting length exactly corresponds to that
            //  of newState.
            //  These distinctions are important if both geometric states do not 
            //  have the same volume.
            virtual pointField transformPoints
            (
                const pointField& undisplacedPoints,
                const deformableBodyGeometricState& newState,
                const bool exactDiameter
            ) const;

        // I/O
            virtual Istream& read(Istream&);
            virtual void write(Ostream&) const;


    // Member Operators

        void operator=(const axisymmetricBodyGeometricState&);


    // IOstream Operators

        friend Istream& operator>>(Istream&, axisymmetricBodyGeometricState&);
        friend Ostream& operator<<(Ostream&, const axisymmetricBodyGeometricState&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
