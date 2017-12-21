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

#include "sampleRBCFieldAxisymmetric.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"

#include "OFstream.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    const word sampleRBCFieldAxisymmetric::RBCFilePrefix = "RBC";

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

defineTypeNameAndDebug(sampleRBCFieldAxisymmetric, 0);
addToRunTimeSelectionTable(functionObject, sampleRBCFieldAxisymmetric, dictionary); 


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::scalarListList>
Foam::sampleRBCFieldAxisymmetric::sampleRBCStats() const
{
    scalarList undefList(4, -VGREAT);
    autoPtr<scalarListList> RBCStatsPtr
    (
        new scalarListList(RBCs_.RBCIndices().size(), undefList)
    );

    scalarListList& RBCStats = RBCStatsPtr();

    forAll(RBCs_.RBCIndices(), i)
    {
        const label RBCI = RBCs_.RBCIndices()[i];
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
Foam::sampleRBCFieldAxisymmetric::computeRBCStats(const label RBCI) const
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
    scalarList result(4);
    result[0] = RBCs_[RBCI].center().x();
    result[1] = sum(RBCField*RBCVolume)/sum(RBCVolume);
    result[2] = min(RBCField);
    result[3] = max(RBCField);

    return result;
}


Foam::word 
Foam::sampleRBCFieldAxisymmetric::sampleDirPath() const
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


void Foam::sampleRBCFieldAxisymmetric::prepare()
{
    if (Pstream::master())
    {
        mkDir(sampleDirPath());
    }
}


void Foam::sampleRBCFieldAxisymmetric::createFile(const label RBCI)
{
    if (Pstream::master())
    {
        sampleFilePtrs_.insert(RBCI, new OFstream(sampleDirPath()/RBCFileName(RBCI)));

        if (debug)
        {
            Info<< "open sample RBC stream: " << sampleFilePtrs_[RBCI]->name() << endl;
        }

        unsigned int w = IOstream::defaultPrecision() + 7;
        *sampleFilePtrs_[RBCI]<< '#' << setw(IOstream::defaultPrecision() + 6)
            << "Time"      << setw(w)
            << "x"         << setw(w)
            << "Mean"      << setw(w) 
            << "Min"       << setw(w) 
            << "Max"       << endl;
    }
}


void Foam::sampleRBCFieldAxisymmetric::closeFile(const label RBCI)
{
    if (debug)
    {
        Info<< "close sample RBC stream: " << sampleFilePtrs_[RBCI]->name() << endl;
    }

    delete sampleFilePtrs_[RBCI];
    sampleFilePtrs_.erase(RBCI);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampleRBCFieldAxisymmetric::sampleRBCFieldAxisymmetric
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
    fieldName_(dict.lookup("field")),
    sampleFilePtrs_()
{
    prepare();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampleRBCFieldAxisymmetric::~sampleRBCFieldAxisymmetric()
{
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

word 
Foam::sampleRBCFieldAxisymmetric::RBCFileName(label i)
{
    std::stringstream sstm;
    sstm << RBCFilePrefix << i;
    return sstm.str();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampleRBCFieldAxisymmetric::read(const dictionary& dict)
{
    dict.lookup("field") >> fieldName_;
    return true;
}

bool Foam::sampleRBCFieldAxisymmetric::execute()
{
    scalarListList RBCStats(sampleRBCStats()());

    if (Pstream::master())
    {
        forAll(RBCs_.RBCIndices(), i)
        {
            const label RBCI = RBCs_.RBCIndices()[i];
            const word fileName = RBCFileName(RBCI);

            if (RBCStats[i][0] > -VGREAT)
            {
                if (!sampleFilePtrs_.found(RBCI))
                {
                    createFile(RBCI);
                }

                unsigned int w = IOstream::defaultPrecision() + 7;

                *sampleFilePtrs_[RBCI]<< ' '                     << setw(w-1)
                                   << time_.timeOutputValue() << setw(w)
                                   << RBCStats[i][0]     << setw(w)
                                   << RBCStats[i][1]     << setw(w)
                                   << RBCStats[i][2]     << setw(w)
                                   << RBCStats[i][3]     << endl;
            }
            else
            {
                // close files for RBCs that left the computational domain
                if (sampleFilePtrs_.found(RBCI))
                {
                    closeFile(RBCI);
                }
            }
        }
        // close files for RBCs that were erased from the RBC list
        forAll(sampleFilePtrs_.toc(), i)
        {
            const label RBCI = sampleFilePtrs_.toc()[i];
            if (findIndex(RBCs_.RBCIndices(), RBCI) == -1)
            {
                    closeFile(RBCI);
            }
        }
    }

    return true;
}


bool Foam::sampleRBCFieldAxisymmetric::write()
{
    return true;
}


bool Foam::sampleRBCFieldAxisymmetric::end()
{
    return true;
}


void Foam::sampleRBCFieldAxisymmetric::updateMesh(const mapPolyMesh& map)
{
}

void Foam::sampleRBCFieldAxisymmetric::movePoints(const polyMesh& mesh)
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
