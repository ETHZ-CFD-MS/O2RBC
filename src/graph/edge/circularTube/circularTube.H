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
    Foam::circularTube

Description
    Abstract class for the description of tubes with circular normal cross 
    section and both an inner and outer diameter.

    The ellipse axis size correspond to the outer diameter.

SourceFiles
    circularTube.C

\*---------------------------------------------------------------------------*/

#ifndef circularTube_H
#define circularTube_H

#include "tubeVisitor.H"

#include "ellipseAxes.H"
#include "Pair.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;

// Forward declaration of friend functions and operators
class circularTube;
Ostream& operator<<(Ostream&, const circularTube&);


/*---------------------------------------------------------------------------*\
                         Class circularTube Declaration
\*---------------------------------------------------------------------------*/

class circularTube
{
    // Private data

        //- Tube axes at end points
        ellipseAxes startAxes_;
        ellipseAxes endAxes_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        circularTube(const circularTube&);

        //- Disallow default bitwise assignment
        void operator=(const circularTube&);

public:

    // Constructors

        //- Construct from components
        circularTube(const ellipseAxes&, const ellipseAxes&);

        //- Construct from ellipseAxes pair
        circularTube(const Pair<ellipseAxes>&);


    //- Destructor
    virtual ~circularTube();


    // Member Functions

        // Access

            const ellipseAxes& startAxes() const
            {
                return startAxes_;
            }

            const ellipseAxes& endAxes() const
            {
                return endAxes_;
            }

        // Setters
            
            virtual void setStartAxes(const ellipseAxes&);
            virtual void setEndAxes(const ellipseAxes&);

        // Geometric queries

            //- Mean outer tube diameter
            virtual scalar meanOuterDiameter() const = 0;

            //- Mean inner tube diameter
            virtual scalar meanInnerDiameter() const = 0;

            //- Mean outer tube radius
            scalar meanOuterRadius() const;

            //- Mean inner tube radius
            scalar meanInnerRadius() const;

            //- Return the tube outer diameter at a given coordinate.
            virtual scalar outerDiameter(const scalar s) const = 0;
    
            //- Return the tube inner diameter at a given coordinate.
            virtual scalar innerDiameter(const scalar s) const = 0;
    
            //- Return the tube outer radius at a given coordinate.
            scalar outerRadius(const scalar s) const;
    
            //- Return the tube inner radius at a given coordinate.
            scalar innerRadius(const scalar s) const;
    
            //- Return the closest ellipseAxes with coordinate < s
            virtual ellipseAxes ellipseAxesBefore(const scalar s) const = 0;

            //- Return the closest ellipseAxes with coordinate >= s
            virtual ellipseAxes ellipseAxesAfter(const scalar s) const = 0;

            //- Return a skeleton of the point tubes
            virtual pointField skeletonPoints() const = 0;

        // Visitor call
            
            //- Call the visitor's visit method (see Visitor pattern).
            virtual void accept(tubeVisitor&) const = 0;

        // Write

            virtual void write(Ostream&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const circularTube&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
