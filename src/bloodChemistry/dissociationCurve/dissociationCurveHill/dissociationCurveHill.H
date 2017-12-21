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
    Foam::dissociationCurveHill

Description
    Dissociation curve based on Hill's equation.

SourceFiles
    dissociationCurveHill.C

\*---------------------------------------------------------------------------*/

#ifndef dissociationCurveHill_H
#define dissociationCurveHill_H

#include "fvCFD.H"

#include "dissociationCurve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class dissociationCurveHill;
Ostream& operator<<(Ostream&, const dissociationCurveHill&);


/*---------------------------------------------------------------------------*\
                Class dissociationCurveHill Declaration
\*---------------------------------------------------------------------------*/

class dissociationCurveHill
    : public dissociationCurve
{
    // Private data

        scalar P50_;

        scalar exponent_;

public:

    //- Runtime type information
    TypeName("Hill");


    // Constructors
    
        //- Construct from dictionary
        dissociationCurveHill
        (
            const dictionary& dict
        );

        //- Construct from components
        dissociationCurveHill
        (
            const scalar P50,
            const scalar exponent
        );

        //- Copy constructor
        dissociationCurveHill
        (
            const dissociationCurveHill&
        );

    // Destructor
        virtual ~dissociationCurveHill()
        {}

    // Member functions
    
        // Evaluation functions
           
            //- compute equilibrium PO2 from saturation
            virtual scalar equilibriumPO2
            (
                const scalar& S
            ) const;

            virtual dimensionedScalar equilibriumPO2
            (
                const dimensionedScalar& S
            ) const;

            virtual tmp<scalarField> equilibriumPO2
            (
                const scalarField& S
            ) const;

            virtual tmp<volScalarField> equilibriumPO2
            (
                const volScalarField& S
            ) const;

            //- compute equilibrium saturation from PO2
            virtual scalar equilibriumS
            (
                const scalar& P
            ) const;

            virtual dimensionedScalar equilibriumS
            (
                const dimensionedScalar& PO2
            ) const;

            virtual tmp<scalarField> equilibriumS
            (
                const scalarField& PO2
            ) const;

            virtual tmp<volScalarField> equilibriumS
            (
                const volScalarField& PO2
            ) const;

            //- compute reaction term
            virtual tmp<volScalarField> reactionTerm
            (
                const volScalarField& PO2,
                const volScalarField& S,
                const volScalarField& inRBC
            ) const;

            virtual tmp<volScalarField> reactionTerm
            (
                const volScalarField& PO2,
                const volScalarField& S
            ) const;

            //- compute derivative of reaction function with respect to PO2
            virtual tmp<volScalarField> ddP
            (
                const volScalarField& PO2,
                const volScalarField& S,
                const volScalarField& inRBC
            ) const;

            virtual tmp<volScalarField> ddS
            (
                const volScalarField& PO2,
                const volScalarField& S
            ) const;


        //- Write
            void write(Ostream&) const;

    //  Ostream operators
            friend Ostream& operator<<
            (
                Ostream&,
                const dissociationCurveHill&
            );

};

Ostream& operator<<
(
    Ostream&,
    const dissociationCurveHill&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
