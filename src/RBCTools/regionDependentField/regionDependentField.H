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
    Foam::regionDependentField

Description

    Describes a volScalarField that takes constant values in different regions
    described with an instance of vascularGraphRegions.
    Additionally, the volScalarField may take different values within RBCs.

SourceFiles
    regionDependentField.C

\*---------------------------------------------------------------------------*/

#ifndef regionDependentField_H
#define regionDependentField_H

#include "fvCFD.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;
class vascularGraphRegions;

// Forward declaration of friend functions and operators
class regionDependentField;
Ostream& operator<<(Ostream&, const regionDependentField&);


/*---------------------------------------------------------------------------*\
                         Class regionDependentField Declaration
\*---------------------------------------------------------------------------*/

class regionDependentField
{
public:

    // Public data types
    
        //- Enumeration for the field averaging method
        enum averagingMethod
        {
            ARITHMETIC,
            HARMONIC
        };


private:

    // Private data
    
        //- Geometric scalar field object
        volScalarField field_;

        //- Reference to vascular graph regions
        const vascularGraphRegions& regions_;

        scalar RBCValue_;
        scalar plasmaValue_;
        scalar wallValue_;
        scalar tissueValue_;

        //- Averaging method
        averagingMethod averagingMethod_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionDependentField(const regionDependentField&);

        //- Disallow default bitwise assignment
        void operator=(const regionDependentField&);

        word averagingMethodToWord(const averagingMethod&) const;

        averagingMethod wordToAveragingMethod(const word&) const;

        scalar weightedHarmonicMean
        (
            const scalar w1, const scalar x1,
            const scalar w2, const scalar x2,
            const scalar w3, const scalar x3

        ) const
        {
            return (w1 + w2 + w3)/(w1/x1 + w2/x2 + w3/x3);
        }

        scalar weightedHarmonicMean
        (
            const scalar w1, const scalar x1,
            const scalar w2, const scalar x2

        ) const
        {
            return (w1 + w2)/(w1/x1 + w2/x2);
        }


public:

    //- Runtime type information
    TypeName("regionDependentField");

    // Constructors

        //- Construct and read given IOobject
        regionDependentField
        (
            const IOobject&,
            const fvMesh&,
            const dictionary& dict,
            const vascularGraphRegions&
        );

    //- Destructor
    virtual ~regionDependentField();


    // Member Functions

        // Update
        
        // Update without taking RBCs into account
        void update();

        // Update by taking RBCs into account
        void update(const volScalarField& inRBC);

        // Access
        
        volScalarField& field()
        {
            return field_;
        }

        const volScalarField& field() const
        {
            return field_;
        }

        scalar RBCValue() const
        {
            return RBCValue_;
        }

        scalar plasmaValue() const
        {
            return plasmaValue_;
        }

        scalar wallValue() const
        {
            return wallValue_;
        }

        scalar tissueValue() const
        {
            return tissueValue_;
        }

        // Write
            bool writeData(Ostream&) const;

            bool write() const;

    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const regionDependentField&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
