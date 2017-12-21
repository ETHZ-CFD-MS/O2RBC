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
    Foam::RBCPath

Description
    Stores information of a RBC path through a graph.

SourceFiles
    RBCPath.C

\*---------------------------------------------------------------------------*/

#ifndef RBCPath_H
#define RBCPath_H

#include "scalarList.H"
#include "labelList.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;
class dictionary;

// Forward declaration of friend functions and operators
class RBCPath;
Ostream& operator<<(Ostream&, const RBCPath&);


/*---------------------------------------------------------------------------*\
                         Class RBCPath Declaration
\*---------------------------------------------------------------------------*/

class RBCPath
{
    // Private data

        //- Path index
        label index_;
        
        //- List of time points
        scalarList pathTimes_;

        //- List of edge labels
        labelList edges_;

        //- List of curvilinear coordinates
        scalarList sCoords_;


    // Private Member Functions
    
        //- Read from files
        void checkLists();

        //- Disallow default bitwise copy construct
        RBCPath(const RBCPath&);

        //- Disallow default bitwise assignment
        void operator=(const RBCPath&);


public:

    // Constructors
        
        //- Construct from dictionary
        RBCPath
        (
            const dictionary&
        );

        //- Construct from components
        RBCPath
        (
            const label index,
            const scalarList& pathTimes,
            const labelList& edges,
            const scalarList& sCoords
        );


    // Member Functions

        // Access

            //- Constant access to index
            label index() const
            {
                return index_;
            }

            //- Constant access to pathTimes
            const scalarList& pathTimes() const
            {
                return pathTimes_;
            }
        
            //- Constant access to edges
            const labelList& edges() const
            {
                return edges_;
            }
        
            //- Constant access to curvilinear coordinate
            const scalarList& sCoords() const
            {
                return sCoords_;
            }

            //- Scale the i-th curvilinear coordinate by a given factor
            void scaleSCoord(const label i, const scalar);
        
            //- Set the i-th edge value
            void setEdge(const label i, const label);
        
            //- Return whether the current path is active, i.e., whether the
            //  current time is within the path time interval.
            bool active(const scalar currentTime) const;

            //- Return whether the current path becomes active at the current
            //  time step
            bool newActive
            (
                const scalar currentTime, 
                const scalar oldTime
            ) const;

            //- Return whether the current path becomes inactive at the current
            //  time step
            bool newInactive
            (
                const scalar currentTime, 
                const scalar oldTime
            ) const;


        // Edit

        // Write

            //- Write to Ostream in simple FOAM format
            void write(Ostream&) const;
            


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const RBCPath&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
