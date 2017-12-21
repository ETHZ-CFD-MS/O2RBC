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

\*---------------------------------------------------------------------------*/

#include "RBCProbe.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"

#include "RBCPath.H"
#include "OFstream.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- comparison operator for RBC statistics
class isNotEqOp
{
public:

    void operator()(scalarList& x, const scalarList& y) const
    {
        const scalar unsetVal(-VGREAT);

        if (x[0] == unsetVal && y[0] != unsetVal)
        {
            x = y;
        }
    }
};

}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(RBCProbe, 0);
addToRunTimeSelectionTable(functionObject, RBCProbe, dictionary); 


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::scalarListList>
Foam::RBCProbe::sampleRBCData() const
{
    autoPtr<scalarListList> RBCDataPtr
    (
        new scalarListList(sampleGraphCoords_.size(), List<scalar>(5, -VGREAT))
    );

    scalarListList& RBCData = RBCDataPtr();

    forAll(sampleGraphCoords_, coordI)
    {
        indexTimeList overlappingRBCsAndTimes =
            probeOverlappingRBCsAndTimesOnProcessor(sampleGraphCoords_[coordI]);

        forAll(overlappingRBCsAndTimes, i)
        {
            label RBCI = overlappingRBCsAndTimes[i].first();
            scalar hittingTime = overlappingRBCsAndTimes[i].second();
            scalarList RBCStats = computeRBCStats(RBCI);
            RBCData[coordI][0] = RBCI;
            RBCData[coordI][1] = hittingTime;
            RBCData[coordI][2] = RBCStats[0];
            RBCData[coordI][3] = RBCStats[1];
            RBCData[coordI][4] = RBCStats[2];
        }
    }
    
    Pstream::listCombineGather(RBCData, isNotEqOp());
    Pstream::listCombineScatter(RBCData);

    return RBCDataPtr;
}
 

Foam::RBCProbe::indexTimeList
Foam::RBCProbe::probeOverlappingRBCsAndTimesOnProcessor
(
    const graphCoordinate& gcSampled
) const
{
    labelList activePathIdx = RBCPaths_.activeIndices();
    indexTimeList overlappingRBCsAndTimes;
    forAll(activePathIdx, i)
    {
        const RBCPath& path = RBCPaths_(activePathIdx[i]);
        if (RBCs_.RBCHasNonEmptyZone(path.index()))
        {
            graphCoordinate gcNow = graphInterpolator_.interpolate
            (
                path.pathTimes(),
                path.edges(),
                path.sCoords(),
                time_.timeOutputValue()
            );

            const scalar previousTime = time_.timeOutputValue() - time_.deltaTValue();

            if (previousTime >= path.pathTimes()[0])
            {
                graphCoordinate gcPrevious = graphInterpolator_.interpolate
                (
                    path.pathTimes(),
                    path.edges(),
                    path.sCoords(),
                    previousTime
                );

                if 
                (
                    graphInterpolator_.isBetweenUpToNeighbourEdge
                    (
                        gcSampled, 
                        gcNow, 
                        gcPrevious
                    )
                )
                {
                    scalar hittingTime = graphInterpolator_.computeGraphCoordinatePassingTime
                                         (
                                             path.pathTimes(),
                                             path.edges(),
                                             path.sCoords(),
                                             gcSampled
                                         );
                    overlappingRBCsAndTimes.append
                    (
                        Tuple2<label, scalar>(path.index(), hittingTime)
                    );
                }
            }
        }
    }

    return overlappingRBCsAndTimes;
}


Foam::scalarList
Foam::RBCProbe::computeRBCStats(const label RBCI) const
{
    const volScalarField& f = RBCs_.RBCMesh().lookupObject<volScalarField>(fieldName_);
    const scalarField& RBCField = RBCs_.RBCMesh().getFieldOnZone
        (
            f, 
            RBCs_.RBCZoneID(RBCI)
        );
    const scalarField& RBCVolume = RBCs_.RBCMesh().getFieldOnZone
        (
            RBCs_.RBCMesh().V().field(), 
            RBCs_.RBCZoneID(RBCI)
        );
    scalarList result(3);
    result[0] = sum(RBCField*RBCVolume)/sum(RBCVolume);
    result[1] = min(RBCField);
    result[2] = max(RBCField);

    return result;
}


void Foam::RBCProbe::prepare()
{
    if (Pstream::master())
    {
        fileName sampleDir;
        fileName sampleSubDir = functionObject::name();

        sampleSubDir = "postProcessing"/sampleSubDir/time_.timeName();

        if (Pstream::parRun())
        {
            sampleDir = time_.path()/".."/sampleSubDir;
        }
        else
        {
            sampleDir = time_.path()/sampleSubDir;
        }

        mkDir(sampleDir);
        sampleFilePtr_ = new OFstream(sampleDir/fieldName_);

        unsigned int w = IOstream::defaultPrecision() + 7;
        *sampleFilePtr_<< '#' << setw(IOstream::defaultPrecision() + 6)
            << "Time" << setw(w)
            << "Edge index" << setw(w)
            << "sCoord" << setw(w)
            << "RBC index" << setw(w)
            << "Hit time" << setw(w)
            << "Mean" << setw(w) 
            << "Min" << setw(w) 
            << "Max" << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBCProbe::RBCProbe
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
:
    functionObject(name),
    time_(time),
    RBCs_
    (
        time_.lookupObject<RBCCollection>("RBCCollection")
    ),
    RBCPaths_
    (
        time_.lookupObject<RBCPathCollection>("RBCPaths")
    ),
    graph_
    (
        time_.lookupObject<geometricEdgeGraph>
        (
            dict.lookup("graphName")
        )
    ),
    graphInterpolator_(graph_),
    sampleGraphCoords_(dict.lookup("sampleGraphCoords")),
    fieldName_(dict.lookup("field")),
    sampleFilePtr_()
{
    prepare();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBCProbe::~RBCProbe()
{
    delete sampleFilePtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RBCProbe::read(const dictionary& dict)
{
    dict.lookup("sampleGraphCoords") >> sampleGraphCoords_;
    dict.lookup("field") >> fieldName_;

    return true;
}

bool Foam::RBCProbe::execute()
{
    scalarListList RBCData(sampleRBCData()());

    if (Pstream::master())
    {
        forAll(sampleGraphCoords_, coordI)
        {
            if (RBCData[coordI][0] > -VGREAT)
            {
                unsigned int w = IOstream::defaultPrecision() + 7;

                *sampleFilePtr_<< ' '                                 << setw(w-1)
                               << time_.timeOutputValue()             << setw(w)
                               << sampleGraphCoords_[coordI].edgeIndex() << setw(w)
                               << sampleGraphCoords_[coordI].sCoord() << setw(w)
                               << RBCData[coordI][0]                  << setw(w)
                               << RBCData[coordI][1]                  << setw(w)
                               << RBCData[coordI][2]                  << setw(w)
                               << RBCData[coordI][3]                  << setw(w)
                               << RBCData[coordI][4]                  << endl;

            }
        }
    }

    return true;
}


bool Foam::RBCProbe::write()
{
    return true;
}


bool Foam::RBCProbe::end()
{
    return true;
}


void Foam::RBCProbe::updateMesh(const mapPolyMesh& map)
{
    // execute(false);
}

void Foam::RBCProbe::movePoints(const polyMesh& mesh)
{
    // execute(false);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
