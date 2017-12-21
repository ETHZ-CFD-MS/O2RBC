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
    Foam::diameterFunctionPries2005

Description
    Computes the RBC diameter from the vessel diameter based on Pries'
    function for the effective ESL width.

SourceFiles
    diameterFunctionPries2005.C

\*---------------------------------------------------------------------------*/

#ifndef diameterFunctionPries2005_H
#define diameterFunctionPries2005_H

#include "fvCFD.H"

#include "diameterFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class diameterFunctionPries2005 Declaration
\*---------------------------------------------------------------------------*/

class diameterFunctionPries2005
    : public diameterFunction
{
    // Private data

        //- Diameter under which the ESL width is zero
        scalar dOff_;

        //- Critical diameter for the three-shape function
        scalar dCrit_;

        //- Top diameter for the three-shape function
        scalar dTop_;

        //-
        scalar d50_;

        scalar eAmp_;

        scalar eWidth_;

        scalar ePeak_;

        //- Maximal ESL width
        scalar wMax_;

    // Private member functions

        //- Compute the effective ESL thickness
        scalar ESLThickness(const scalar diameter) const;


public:

    //- Runtime type information
    TypeName("Pries2005");


    // Constructors
    
        // Construct from dictionary
        diameterFunctionPries2005
        (
            const dictionary& dict
        );

    // Destructor
        virtual ~diameterFunctionPries2005()
        {}

    // Member functions
    
        // Evaluation functions
           
            //- Evaluates the RBC diameter from the vessel diameter
            scalar operator() (const scalar) const;

        //- Write
            void write(Ostream&) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
