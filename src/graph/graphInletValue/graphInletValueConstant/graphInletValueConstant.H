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
    Foam::graphInletValueConstant

Description
    Returns a constant value at an inlet edge of a geometric graph.

SourceFiles
    graphInletValueConstant.C

\*---------------------------------------------------------------------------*/

#ifndef graphInletValueConstant_H
#define graphInletValueConstant_H

#include "fvCFD.H"

#include "graphInletValue.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class geometricEdgeGraph;

/*---------------------------------------------------------------------------*\
                Class graphInletValueConstant Declaration
\*---------------------------------------------------------------------------*/

class graphInletValueConstant
    : public graphInletValue
{
    // Private data

        //- Map between edge indices and constant values
        Map<scalar> inletValues_;

public:

    //- Runtime type information
    TypeName("constant");


    // Constructors
    
        // Construct from dictionary
        graphInletValueConstant
        (
            const dictionary& dict,
            const geometricEdgeGraph& graph
        );

    // Destructor
        virtual ~graphInletValueConstant()
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


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
