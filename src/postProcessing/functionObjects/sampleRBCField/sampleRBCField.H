/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Description
    Function object for sampling of fields in RBCs that are currently in the
    computational domain.

    Typical entry in controlDict:

    objectName
    {
        type        sampleRBCField;
        field       Hb;
    }

SourceFiles
    sampleRBCField.C

\*---------------------------------------------------------------------------*/

#ifndef sampleRBCField_H
#define sampleRBCField_H

#include "functionObject.H"

#include <fstream>

#include "RBCCollection.H"
#include "RBCPathCollection.H"
#include "geometricEdgeGraph.H"
#include "graphCoordinateInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    bool fileExists(const char *fileName);

/*---------------------------------------------------------------------------*\
                      Class sampleRBCField Declaration
\*---------------------------------------------------------------------------*/

class sampleRBCField
    :
        public functionObject
{
    // Private data

        //- Time
        const Time& time_; 

        //- RBC collection
        const RBCCollection& RBCs_;

        //- RBC paths
        const RBCPathCollection& RBCPaths_;

        //- Name of the field to sample
        word fieldName_;

        //- Graph
        const geometricEdgeGraph& graph_;

        //- Interpolation object for the RBC paths
        graphCoordinateInterpolation graphInterpolator_;

        //- Current open files
        List<std::ofstream*> sampleFilePtrs_;

    // Private member functions

        //- Return a list of RBC statistics for each RBC path
        virtual autoPtr<scalarListList> sampleRBCStats() const;

        //- Compute the RBC statistics for a given processor
        virtual scalarList computeRBCStats(const label RBCI) const;

        //- Create the output directory path (uses the startTime)
        word sampleDirPath() const;

        //- Prepare output file and its header
        void prepare();

        //- Create a file for a given index
        void openFile(const label, const bool writeHeader);

        //- Close a file for a given index
        void closeFile(const label);
        

public:

    //- Runtime type information
    TypeName("sampleRBCField");

    //- Static data

        //- Prefix for RBC file
        static const word RBCFilePrefix;


    // Static member functions
        
        // Return the RBC file name
        static word RBCFileName(label);


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        sampleRBCField
        (
            const word& name,
            const Time& time,
            const dictionary& dict
        );

    //- Destructor
     ~sampleRBCField();


    // Member Functions

        //- Access time
        const Time& getTime() const
        {
            return time_; 
        }

        //- Read the system calls
         bool read(const dictionary&);

        //- Called at each ++ or += of the time-loop. 
         bool execute();

        //- Called at each ++ or += of the time-loop. 
         bool write();

        //- Called when Time::run() determines that the time-loop exits
         bool end();

        //- Update for changes of mesh
         void updateMesh(const mapPolyMesh&);

        //- Update for changes of mesh
         void movePoints(const polyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
