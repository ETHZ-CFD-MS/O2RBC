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

#include <iomanip>

#include "sampleRBCField.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"

#include "RBCPath.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    bool fileExists(const char *fileName)
    {
        std::ifstream infile(fileName);
        return infile.good();
    }

    const word sampleRBCField::RBCFilePrefix = "RBC";

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

defineTypeNameAndDebug(sampleRBCField, 0);
addToRunTimeSelectionTable(functionObject, sampleRBCField, dictionary); 


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::scalarListList>
Foam::sampleRBCField::sampleRBCStats() const
{
    scalarList undefList(5, -VGREAT);
    autoPtr<scalarListList> RBCStatsPtr
    (
        new scalarListList(RBCPaths_.size(), undefList)
    );

    scalarListList& RBCStats = RBCStatsPtr();

    forAll(RBCPaths_, i)
    {
        const RBCPath& cPath = RBCPaths_[i];
        const label RBCI = cPath.index();
        if (RBCs_.RBCHasNonEmptyZone(RBCI))
        {
            RBCStats[i] = computeRBCStats(RBCI);
        }
        else
        {
            RBCStats[i] = undefList;
        }
    }
    
    Pstream::listCombineGather(RBCStats, isNotEqOp());
    Pstream::listCombineScatter(RBCStats);

    return RBCStatsPtr;
}
 

Foam::scalarList
Foam::sampleRBCField::computeRBCStats(const label RBCI) const
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
    const RBCPath& path = RBCPaths_(RBCI);
    graphCoordinate gc = graphInterpolator_.interpolate
    (
        path.pathTimes(),
        path.edges(),
        path.sCoords(),
        time_.timeOutputValue()
    );

    scalarList result(5);
    result[0] = sum(RBCField*RBCVolume)/sum(RBCVolume);
    result[1] = min(RBCField);
    result[2] = max(RBCField);
    result[3] = gc.edgeIndex();
    result[4] = gc.sCoord();

    return result;
}


Foam::word 
Foam::sampleRBCField::sampleDirPath() const
{
    fileName sampleDir;
    fileName sampleSubDir = functionObject::name();
    sampleSubDir = "postProcessing"/sampleSubDir/
                    time_.timeName(time_.timeToUserTime(time_.startTime().value()));

    if (Pstream::parRun())
    {
        sampleDir = time_.path()/".."/sampleSubDir;
    }
    else
    {
        sampleDir = time_.path()/sampleSubDir;
    }
    return sampleDir;
}

void Foam::sampleRBCField::prepare()
{
    if (Pstream::master())
    {
        mkDir(sampleDirPath());
    }
}


void Foam::sampleRBCField::openFile(const label i, const bool writeHeader)
{
    if (Pstream::master())
    {
        const label RBCPathI = RBCPaths_[i].index();
        const fileName RBCFilePath = sampleDirPath()/RBCFileName(RBCPathI);
        sampleFilePtrs_[i] = new std::ofstream(RBCFilePath.c_str(), std::ios_base::app);

        if (debug)
        {
            Info<< "open sample RBC stream: " << RBCFilePath << endl;
        }

        if (writeHeader)
        {
            unsigned int w = IOstream::defaultPrecision() + 7;
            *sampleFilePtrs_[i]<< '#' << std::setw(IOstream::defaultPrecision() + 6)
                << "Time"      << std::setw(w)
                << "Mean"      << std::setw(w) 
                << "Min"       << std::setw(w) 
                << "Max"       << std::setw(w) 
                << "EdgeIndex" << std::setw(w) 
                << "sCoord"    << std::endl;
        }
    }
}


void Foam::sampleRBCField::closeFile(const label i)
{
    if (debug)
    {
        const label RBCPathI = RBCPaths_[i].index();
        const fileName RBCFilePath = sampleDirPath()/RBCFileName(RBCPathI);
        Info<< "close sample RBC stream: " << RBCFilePath << endl;
    }

    delete sampleFilePtrs_[i];
    sampleFilePtrs_[i] = 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampleRBCField::sampleRBCField
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
    fieldName_(dict.lookup("field")),
    graph_
    (
        time_.lookupObject<geometricEdgeGraph>
        (
            dict.lookup("graphName")
        )
    ),
    graphInterpolator_(graph_),
    sampleFilePtrs_(RBCPaths_.size(), NULL)
{
    prepare();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampleRBCField::~sampleRBCField()
{
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

word 
Foam::sampleRBCField::RBCFileName(label i)
{
    std::stringstream sstm;
    sstm << RBCFilePrefix << i;
    return sstm.str();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampleRBCField::read(const dictionary& dict)
{
    dict.lookup("field") >> fieldName_;

    return true;
}

bool Foam::sampleRBCField::execute()
{
    scalarListList RBCStats(sampleRBCStats()());

    if (Pstream::master())
    {
        forAll(RBCPaths_, i)
        {
            const label RBCPathI = RBCPaths_[i].index();
            const fileName RBCFilePath = sampleDirPath()/RBCFileName(RBCPathI);

            if (RBCStats[i][0] > -VGREAT)
            {
                if (!sampleFilePtrs_[i])
                {
                    bool writeHeader = true;
                    if (fileExists(RBCFilePath.c_str()))
                    {
                        writeHeader = false;
                    }
                    openFile(i, writeHeader);
                }

                unsigned int w = IOstream::defaultPrecision() + 7;

                if (!sampleFilePtrs_[i]->is_open())
                {
                    sampleFilePtrs_[i]->open(RBCFilePath, ios_base::app);
                }

                *sampleFilePtrs_[i]<< ' '                     << std::setw(w-1)
                                   << time_.timeOutputValue() << std::setw(w)
                                   << RBCStats[i][0]     << std::setw(w)
                                   << RBCStats[i][1]     << std::setw(w)
                                   << RBCStats[i][2]     << std::setw(w)
                                   << RBCStats[i][3]     << std::setw(w)
                                   << RBCStats[i][4]     << std::endl;
                sampleFilePtrs_[i]->close();
                sampleFilePtrs_[i]->clear();
            }
            else
            {
                if (sampleFilePtrs_[i])
                {
                    closeFile(i);
                }
            }
        }
    }

    return true;
}


bool Foam::sampleRBCField::write()
{
    return true;
}


bool Foam::sampleRBCField::end()
{
    return true;
}


void Foam::sampleRBCField::updateMesh(const mapPolyMesh& map)
{
}

void Foam::sampleRBCField::movePoints(const polyMesh& mesh)
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
