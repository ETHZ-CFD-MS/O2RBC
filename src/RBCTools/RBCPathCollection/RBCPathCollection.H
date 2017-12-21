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
    Foam::RBCPathCollection

Description
    Container for RBCPaths with IO functionality.

SourceFiles
    RBCPathCollection.C

\*---------------------------------------------------------------------------*/

#ifndef RBCPathCollection_H
#define RBCPathCollection_H

#include "IOdictionary.H"
#include "PtrList.H"
#include "labelList.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class RBCPath;


// Forward declaration of friend functions and operators
class RBCPathCollection;
Ostream& operator<<(Ostream&, const RBCPathCollection&);


/*---------------------------------------------------------------------------*\
                         Class RBCPathCollection Declaration
\*---------------------------------------------------------------------------*/

class RBCPathCollection
:
    public IOdictionary  // derives from IOdictionary so that a RBC collection
                         // can be read as a dictionary
{
    // Private data

        //- Map of RBCPath objects
        Map<autoPtr<RBCPath> > RBCPaths_;

        //- List of RBC indices
        labelList RBCIdxList_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        RBCPathCollection(const RBCPathCollection&);

        //- Disallow default bitwise assignment
        void operator=(const RBCPathCollection&);


public:

    // Runtime type information
    ClassName("RBCPathCollection");

    // Constructors

        //- Construct from IOobject by reading RBC paths from file
        RBCPathCollection(const IOobject&);

        //- Construct from IOobject and RBC paths
        RBCPathCollection(const IOobject&, const Map<autoPtr<RBCPath> >&);


    //- Destructor
    ~RBCPathCollection();


    // Member Functions

        // Access

            //- Size of the path collection        
            label size() const;

            //- Return indices of active paths
            labelList activeIndices() const;
            
            //- Return indices of inactive paths
            labelList inactiveIndices() const;
            
            //- Return indices of new active paths, i.e.
            //  those that became active during the last time step
            labelList newActiveIndices() const;
            
            //- Return indices of inactive paths
            //  those that became inactive during the last time step
            labelList newInactiveIndices() const;

            //- Return whether a path with the given index exists
            bool found(const label RBCI) const
            {
                return RBCPaths_.found(RBCI);
            }


        // Write
            bool writeData(Ostream&) const;

    // Member operators
        
        //- Constant access to individual RBC paths using the RBC index
        const RBCPath& operator()(const label RBCI) const;

        //- Non-constant to access individual RBC paths using the RBC index
        RBCPath& operator()(const label RBCI);

        //- Constant access to individual RBC paths using an index from 0 to
        //  size() - 1
        const RBCPath& operator[](const label) const;

        //- Non-constant to access individual RBC paths using an index from 0 to
        //  size() - 1
        RBCPath& operator[](const label);

    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const RBCPathCollection&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

