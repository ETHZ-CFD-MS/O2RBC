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
    Foam::graphInletValue

Description
    Abstract class for an inlet value defined on a graph.

SourceFiles
    graphInletValue.C

\*---------------------------------------------------------------------------*/

#ifndef graphInletValue_H
#define graphInletValue_H

#include "fvCFD.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class graphInletValue;
class geometricEdgeGraph;
    
Ostream& operator<<
(
    Ostream&,
    const graphInletValue&
);

/*---------------------------------------------------------------------------*\
                Class graphInletValue Declaration
\*---------------------------------------------------------------------------*/

class graphInletValue
{

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        graphInletValue(const graphInletValue&);

        //- Disallow default bitwise assignment
        void operator=(const graphInletValue&);

protected:

    // Protected data

        //- Reference to graph
        const geometricEdgeGraph& graph_;

        //- Name of the field defined at the inlet
        word fieldName_;

        //- List of edges assigned with an inlet value
        labelList assignedEdges_;

    // Protected Member Functions

        //- Check whether the leaf edges have an assigned inlet value
        void checkLeafEdges() const;


public:
    
    //- Runtime type information
    TypeName("graphInletValue");
        

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            graphInletValue,
            dictionary,
            (
                const dictionary& dict,
                const geometricEdgeGraph& graph
            ),
            (dict, graph)
        );

    // Selectors

        //- Return a reference to the selected interpolation scheme
        static autoPtr<graphInletValue> New
        (
            const dictionary& dict,
            const geometricEdgeGraph& graph
        );


    // Constructors
    
        // Construct from dictionary
        graphInletValue
        (
            const dictionary& dict,
            const geometricEdgeGraph& graph
        );

    // Destructor
        virtual ~graphInletValue()
        {};

    // Member functions
    
        // Access

            //- Return the field name
            const word& fieldName() const
            {
                return fieldName_;
            }

        // Evaluation operator
           
            //- Return the inlet value at a given edge
            virtual scalar inletValue
            (
                const label edgeIndex,
                const scalar time
            ) const = 0;


        //- Write
            virtual void write(Ostream&) const = 0;

    // IOstream operator
        
        friend Ostream& operator<<(Ostream&, const graphInletValue&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
