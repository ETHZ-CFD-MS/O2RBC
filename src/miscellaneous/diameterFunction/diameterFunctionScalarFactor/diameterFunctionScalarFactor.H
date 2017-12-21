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
    Foam::diameterFunctionScalarFactor

Description
    Computes the RBC diameter from the vessel diameter using a scalar factor.

SourceFiles
    diameterFunctionScalarFactor.C

\*---------------------------------------------------------------------------*/

#ifndef diameterFunctionScalarFactor_H
#define diameterFunctionScalarFactor_H

#include "fvCFD.H"

#include "diameterFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class diameterFunctionScalarFactor Declaration
\*---------------------------------------------------------------------------*/

class diameterFunctionScalarFactor
    : public diameterFunction
{
    // Private data

        //- Factor that yields the RBC diameter from the lumen diameter
        scalar factor_;

        //- Lower bound for the RBC diameter
        scalar minDiameter_;

        //- Upper bound for the RBC diameter
        scalar maxDiameter_;


public:

    //- Runtime type information
    TypeName("scalarFactor");


    // Constructors
    
        // Construct from dictionary
        diameterFunctionScalarFactor
        (
            const dictionary& dict
        );

        // Construct from components
        diameterFunctionScalarFactor
        (
            const scalar factor,
            const scalar minDiameter,
            const scalar maxDiameter
        );


    // Destructor
        virtual ~diameterFunctionScalarFactor()
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
