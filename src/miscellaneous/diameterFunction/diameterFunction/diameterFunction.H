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
    Foam::diameterFunction

Description
    Abstract class for a diameter dependent function.

SourceFiles
    diameterFunction.C

\*---------------------------------------------------------------------------*/

#ifndef diameterFunction_H
#define diameterFunction_H

#include "fvCFD.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class diameterFunction;
    
Ostream& operator<<
(
    Ostream&,
    const diameterFunction&
);

/*---------------------------------------------------------------------------*\
                Class diameterFunction Declaration
\*---------------------------------------------------------------------------*/

class diameterFunction
{

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        diameterFunction(const diameterFunction&);

        //- Disallow default bitwise assignment
        void operator=(const diameterFunction&);


public:
    
    //- Runtime type information
    TypeName("diameterFunction");
        

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            diameterFunction,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );

    // Selectors

        //- Return a reference to the selected interpolation scheme
        static autoPtr<diameterFunction> New
        (
            const dictionary& dict
        );


    // Constructors
    
        // Construct from dictionary
        diameterFunction
        (
            const dictionary& dict
        );

    // Destructor
        virtual ~diameterFunction()
        {};

    // Member functions
    
        // Evaluation operator
           
            //- Evaluates the RBC diameter from the capillary diameter
            virtual scalar operator() (const scalar) const = 0;


        //- Write
            virtual void write(Ostream&) const = 0;

    // IOstream operator
        
        friend Ostream& operator<<(Ostream&, const diameterFunction&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
