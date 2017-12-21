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
    Foam::RBCProperty

Description
    A property of RBCs with templated type.

SourceFiles
    RBCProperty.C

\*---------------------------------------------------------------------------*/

#ifndef RBCProperty_H
#define RBCProperty_H

#include "fvCFD.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class T>
class RBCProperty;

template<class T>
Ostream& operator<<
(
    Ostream&,
    const RBCProperty<T>&
);

/*---------------------------------------------------------------------------*\
                Class RBCProperty Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class RBCProperty
{

private:

    // Private Member Functions

        //- Check the type of property
        void checkType() const;

        //- Disallow default bitwise copy construct
        RBCProperty(const RBCProperty<T>&);

        //- Disallow default bitwise assignment
        void operator=(const RBCProperty<T>&);

protected:

    // Protected data

        //- Type of RBCProperty
        const word type_;

        //- Constant value
        T constantValue_;

        //- Values for edge-defined property
        Map<T> edgeValues_;


public:
    
    //- Runtime type information
    ClassName("RBCProperty");
        
    // Constructors
    
        // Construct from dictionary
        RBCProperty
        (
            const dictionary& dict
        );

    // Destructor
        virtual ~RBCProperty()
        {};

    // Member functions
    
        // Access

            //- Return a list with all stored values
            List<T> values() const;

        // Evaluation operator
           
            //- Return the inlet value at a given edge
            T value(const label RBCIdx, const label edgeI) const; 

        //- Write
            virtual void write(Ostream&) const;

    // IOstream operator
        friend Ostream& operator<< <T>(Ostream&, const RBCProperty<T>&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "RBCProperty.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
