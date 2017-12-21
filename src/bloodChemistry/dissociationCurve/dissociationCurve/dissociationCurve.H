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
    Foam::dissociationCurve

Description
    Foam::dissociationCurve

SourceFiles
    dissociationCurve.C

\*---------------------------------------------------------------------------*/

#ifndef dissociationCurve_H
#define dissociationCurve_H

#include "fvCFD.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class dissociationCurve;
    
Ostream& operator<<
(
    Ostream&,
    const dissociationCurve&
);

/*---------------------------------------------------------------------------*\
                Class dissociationCurve Declaration
\*---------------------------------------------------------------------------*/

class dissociationCurve
{

private:

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const dissociationCurve&);


public:
    
    //- Runtime type information
    TypeName("dissociationCurve");
        

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            dissociationCurve,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );

    // Selectors

        //- Return a reference to the selected dissociation curve
        static autoPtr<dissociationCurve> New
        (
            const dictionary& dict
        );


    // Constructors
    
        // Construct null
        dissociationCurve() {}

        // Construct from dictionary
        dissociationCurve
        (
            const dictionary& dict
        );

        // Copy constructor
        dissociationCurve
        (
            const dissociationCurve&
        );

    // Destructor
        virtual ~dissociationCurve()
        {};

    // Member functions
    
        // Evaluation functions
           
            //- compute equilibrium PO2 from saturation
            virtual scalar equilibriumPO2
            (
                const scalar& S
            ) const = 0;

            virtual dimensionedScalar equilibriumPO2
            (
                const dimensionedScalar& S
            ) const = 0;

            virtual tmp<scalarField> equilibriumPO2
            (
                const scalarField& S
            ) const = 0;

            virtual tmp<volScalarField> equilibriumPO2
            (
                const volScalarField& S
            ) const = 0;

            //- compute equilibrium saturation from PO2
            virtual scalar equilibriumS
            (
                const scalar& P
            ) const = 0;

            virtual dimensionedScalar equilibriumS
            (
                const dimensionedScalar& PO2
            ) const = 0;

            virtual tmp<scalarField> equilibriumS
            (
                const scalarField& PO2
            ) const = 0;

            virtual tmp<volScalarField> equilibriumS
            (
                const volScalarField& PO2
            ) const = 0;

            //- compute reaction term
            virtual tmp<volScalarField> reactionTerm
            (
                const volScalarField& PO2,
                const volScalarField& S,
                const volScalarField& inRBC
            ) const = 0;

            virtual tmp<volScalarField> reactionTerm
            (
                const volScalarField& PO2,
                const volScalarField& S
            ) const = 0;

            //- compute derivative of reaction function with respect to PO2
            virtual tmp<volScalarField> ddP
            (
                const volScalarField& PO2,
                const volScalarField& S,
                const volScalarField& inRBC
            ) const = 0;

            //- compute derivative of reaction function with respect to S
            virtual tmp<volScalarField> ddS
            (
                const volScalarField& PO2,
                const volScalarField& S
            ) const = 0;


        //- Write
            virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
