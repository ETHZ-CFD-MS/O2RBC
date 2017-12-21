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
    Foam::graphInletValueRandom

Description
    Returns a random value at an inlet edge of a geometric graph.

SourceFiles
    graphInletValueRandom.C

\*---------------------------------------------------------------------------*/

#ifndef graphInletValueRandom_H
#define graphInletValueRandom_H

#include "fvCFD.H"

#include "graphInletValue.H"
#include "randomVariable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class geometricEdgeGraph;

/*---------------------------------------------------------------------------*\
                Class graphInletValueRandom Declaration
\*---------------------------------------------------------------------------*/

class graphInletValueRandom
    : public graphInletValue
{
    // Private data

        //- Map between edge indices and random variables
        Map<randomVariable> inletRandomVars_;

public:

    //- Runtime type information
    TypeName("random");


    // Constructors
    
        // Construct from dictionary
        graphInletValueRandom
        (
            const dictionary& dict,
            const geometricEdgeGraph& graph
        );

    // Destructor
        virtual ~graphInletValueRandom()
        {}

    // Member functions
    
        // Evaluation functions
           
            //- Return the inlet value at a given edge
            virtual scalar inletValue
            (
                const label edgeIndex,
                const scalar time
            ) const;

        //- Write
            void write(Ostream&) const;

    // IOstream Operators

        friend Istream& operator>>(Istream&, graphInletValueRandom&);
        friend Ostream& operator<<(Ostream&, const graphInletValueRandom&);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
