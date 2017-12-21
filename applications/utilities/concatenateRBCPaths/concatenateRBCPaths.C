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

Application
    concatenateRBCPaths

Description
    Concatenate RBC paths that are in different files and write them to a 
    single file.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ListOps.H"

#include "RBCPathCollection.H"
#include "RBCPath.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


autoPtr<RBCPath>
mergeRBCPath
(
    const label pathI, 
    const PtrList<RBCPathCollection>& pathCollections
)
{
    labelList pathInCollectionIds;
    forAll(pathCollections, cI)
    {
        if (pathCollections[cI].found(pathI))
        {
            pathInCollectionIds.append(cI);
        }
    }
    label newListSize = 0;
    forAll(pathInCollectionIds, i)
    {
        label cI = pathInCollectionIds[i];
        newListSize += pathCollections[cI](pathI).pathTimes().size();
    }
    scalarList newTimes(newListSize);
    labelList newEdges(newListSize);
    scalarList newSCoords(newListSize);
    label offset = 0;
    forAll(pathInCollectionIds, i)
    {
        label cI = pathInCollectionIds[i];
        const RBCPath& cPath = pathCollections[cI](pathI);
        forAll(cPath.pathTimes(), j)
        {
            newTimes[offset + j] = cPath.pathTimes()[j];
            newEdges[offset + j] = cPath.edges()[j];
            newSCoords[offset + j] = cPath.sCoords()[j];
        }
        offset += cPath.pathTimes().size();
    }
    return autoPtr<RBCPath>
        (
            new RBCPath
            (
                pathI,
                newTimes,
                newEdges,
                newSCoords
            )
        );
}


autoPtr<RBCPathCollection>
mergeRBCPathCollections
(
    const PtrList<RBCPathCollection>& pathCollections,
    const word& outputFile
)
{
    Map<autoPtr<RBCPath> > RBCPaths;
    labelList visitedRBCPathIds;
    forAll(pathCollections, cI)
    {
        forAll(pathCollections[cI], pI)
        {
            const RBCPath& cPath = pathCollections[cI][pI];
            if (findIndex(visitedRBCPathIds, cPath.index()) < 0)
            {
                RBCPaths.insert
                (
                    cPath.index(),
                    mergeRBCPath(cPath.index(), pathCollections)
                );
                visitedRBCPathIds.append(cPath.index());
            }
        }
    }

    const Time& runTime = pathCollections[0].time();

    autoPtr<RBCPathCollection> pathCollectionPtr
    (
        new RBCPathCollection
        (
            IOobject
            (
                outputFile,
                runTime.caseConstant(),
                runTime,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            RBCPaths
        )
    );
    return pathCollectionPtr;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Concatenate multiple RBC path dictionaries and write them "
        "to a single file."
    );
    #include "addRegionOption.H"
    argList::addOption
    (
        "fileNames",
        "file names",
        "Names of the RBC path dictionaries to concatenate"
    );
    argList::addOption
    (
        "outputFile",
        "output file",
        "Name of the output dictionary"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    wordList fileNames(args.optionReadList<word>("fileNames"));
    word outputFile(args.optionRead<word>("outputFile"));

    PtrList<RBCPathCollection> pathCollections;

    forAll(fileNames, i)
    {
        pathCollections.append
        (
            new RBCPathCollection
            (
                IOobject
                (
                    fileNames[i],
                    runTime.constant(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false // do not register in database
                )
            )
        );
    }
    autoPtr<RBCPathCollection> mergedCollectionPtr =
        mergeRBCPathCollections(pathCollections, outputFile);
    mergedCollectionPtr->regIOobject::write();

    Info<< "Wrote merged RBC path collections." << nl
        << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

