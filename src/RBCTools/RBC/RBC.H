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
    Foam::RBC

Description
    Represents a RBC with the associated geometrical information.

SourceFiles
    RBC.C

\*---------------------------------------------------------------------------*/

#ifndef RBC_H
#define RBC_H

#include "fvCFD.H"
#include "deformableBodyGeometricState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class RBC;
Istream& operator>>(Istream&, RBC&);
Ostream& operator<<(Ostream&, const RBC&);

/*---------------------------------------------------------------------------*\
                             Class RBC Declaration
\*---------------------------------------------------------------------------*/


class RBC
{
    // Private data

        //- RBC identifier
        label index_;

        //- RBC geometric data
        autoPtr<deformableBodyGeometricState> geometricStatePtr_;

        //- Name of sample mesh 
        word sampleMeshName_;

        //- Initial PO2
        scalar initialPO2_;

        //- Inlet PO2
        scalar inletPO2_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const RBC&);


public:
    
    // Constructors

        // Construct from components
        RBC
        (
            const label index,
            const deformableBodyGeometricState& geometricState,
            const word& sampleMeshName,
            const scalar initialPO2
        );

        // Construct from dictionary
        RBC
        (
            const dictionary&
        );

        // Copy constructor
        RBC
        (
            const RBC&
        );

    // Member functions

        // Access
        
            //- Access to the RBC index
            label index() const
            {
                return index_;
            }
        
            const deformableBodyGeometricState& geometricState() const
            {
                return geometricStatePtr_();
            }

            const word sampleMeshName() const
            {
                return sampleMeshName_;
            }

            scalar initialPO2() const
            {
                return initialPO2_;
            }

            scalar inletPO2() const
            {
                return inletPO2_;
            }

            //- Access to the RBC center
            point center() const;

        // Edit
        
            void setInletPO2(const scalar inletPO2)
            {
                inletPO2_ = inletPO2;
            }
            
            void setGeometricState(const deformableBodyGeometricState&);

        
        //- Write
            void write(Ostream&) const;

        // IOstream Operators

            friend Istream& operator>>(Istream&, RBC&);
            friend Ostream& operator<<(Ostream&, const RBC&);

};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
