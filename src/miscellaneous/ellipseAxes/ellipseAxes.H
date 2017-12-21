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
    Foam::ellipseAxes

Description
    Holds geometric information for a two-dimensional ellipse placed in 3D space.

SourceFiles
    ellipseAxes.C

\*---------------------------------------------------------------------------*/

#ifndef ellipseAxes_H
#define ellipseAxes_H

#include "scalar.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;
class ellipseAxes;

// Forward declaration of friend functions and operators
Ostream& operator<<(Ostream&, const ellipseAxes&);

/*---------------------------------------------------------------------------*\
                         Class ellipseAxes Declaration
\*---------------------------------------------------------------------------*/

class ellipseAxes
{
    // Private data

        //- Ellipse axes
        vector minorAxis_;
        vector majorAxis_;

        //- Ellipse radii
        scalar minorRadius_;
        scalar majorRadius_;

public:

    // Constructors

        // Construct null
        ellipseAxes();

        //- Construct from components
        ellipseAxes
        (
            const vector& minorAxis,
            const vector& majorAxis,
            const scalar minorRadius,
            const scalar majorRadius
        );


        //- Construct as copy
        ellipseAxes(const ellipseAxes&);


    //- Destructor
    ~ellipseAxes();


    // Member Functions

        // Access
            const vector& minorAxis() const
            {
                return minorAxis_;
            }

            const vector& majorAxis() const
            {
                return majorAxis_;
            }

            scalar minorRadius() const
            {
                return minorRadius_;
            }

            scalar majorRadius() const
            {
                return majorRadius_;
            }

        // Edit

            void setAxes
            (
                const vector& minorAxis, 
                const vector& majorAxis
            );

            void setRadii
            (
                const scalar minorRadius,
                const scalar majorRadius
            );

        // Interpolation

            //- Interpolate between current object and the axes passed as
            //  argument. The interpolation factor lambda should be between 0
            //  and 1.
            //  This feature is still experimental...
            ellipseAxes interpolate
            (
                const scalar lambda, 
                const ellipseAxes&
            ) const;


        // Check
            void checkAxesAndRadii() const;



    // Member Operators

        void operator=(const ellipseAxes&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const ellipseAxes&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
