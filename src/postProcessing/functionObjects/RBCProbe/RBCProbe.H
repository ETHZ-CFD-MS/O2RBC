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
    Function object for sampling of fields on RBCs that pass over a given 
    position in a geometric graph.

    Typical entry in controlDict:

    objectName
    {
        type        RBCProbe;
        graphName   graphDict;
        field       Hb;
        sampleGraphCoords
        ( 
            (0  2e-6) 
            (7 38e-6) 
        );
    }

SourceFiles
    RBCProbe.C

\*---------------------------------------------------------------------------*/

#ifndef RBCProbe_H
#define RBCProbe_H

#include "functionObject.H"

#include "RBCCollection.H"
#include "RBCPathCollection.H"
#include "geometricEdgeGraph.H"
#include "graphCoordinate.H"
#include "graphCoordinateInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         A templated functionObject
\*---------------------------------------------------------------------------*/

class RBCProbe
    :
        public functionObject
{
    // Public typedefs
    public:
        typedef List<Tuple2<label, scalar> > indexTimeList;

    private:
    // Private data

        //- Time
        const Time& time_; 

        //- RBC collection
        const RBCCollection& RBCs_;

        //- RBC paths
        const RBCPathCollection& RBCPaths_;

        //- Graph
        const geometricEdgeGraph& graph_;

        //- Interpolation object for the RBC paths
        graphCoordinateInterpolation graphInterpolator_;

        //- List of graph coordinates to sample
        List<graphCoordinate> sampleGraphCoords_;

        //- Name of the field to sample
        word fieldName_;

        //- File pointer
        OFstream* sampleFilePtr_;

    // Private member functions

        //- Return RBC data for each probe
        virtual autoPtr<scalarListList> sampleRBCData() const;

        //- Return the RBCs that overlapped a probe on a given processor
        virtual indexTimeList
        probeOverlappingRBCsAndTimesOnProcessor
        (
            const graphCoordinate&
        ) const;

        //- Compute the RBC statistics for a given processor
        virtual scalarList computeRBCStats(const label RBCI) const;

        //- Prepare output file and its header
        void prepare();
        

public:

    //- Runtime type information
    TypeName("RBCProbe");

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        RBCProbe
        (
            const word& name,
            const Time& time,
            const dictionary& dict
        );

    //- Destructor
     ~RBCProbe();


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
