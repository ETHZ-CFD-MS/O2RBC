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
    concatenateEdgeVelocities

Description
    Concatenate edge velocities that are in difference files and write them to
    a single file.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ListOps.H"

using namespace Foam;

typedef Tuple2<scalar, SLList<scalar> > timeVelocities;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


List<timeVelocities>
mergeEdgeVelocities
(
    const PtrList<List<timeVelocities> > edgeVelocitiesPtrList
)
{
    label nStepTotal = 0;
    forAll(edgeVelocitiesPtrList, i)
    {
        nStepTotal += edgeVelocitiesPtrList[i].size();
    }
    List<timeVelocities> mergedList(nStepTotal);
    label offset = 0;
    forAll(edgeVelocitiesPtrList, i)
    {
        forAll(edgeVelocitiesPtrList[i], timeI)
        {
            mergedList[timeI + offset] = edgeVelocitiesPtrList[i][timeI];
        }
        offset += edgeVelocitiesPtrList[i].size();
    }
    return mergedList;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Concatenate edge velocities from different files and write them "
        "to a single file."
    );
    #include "addRegionOption.H"
    argList::addOption
    (
        "fileNames",
        "file names",
        "Names of the edge velocity files to concatenate"
    );
    argList::addOption
    (
        "outputFile",
        "output file",
        "Name of the output file"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    wordList fileNames(args.optionReadList<word>("fileNames"));
    word outputFile(args.optionRead<word>("outputFile"));

    PtrList<List<timeVelocities> > edgeVelocitiesPtrList;

    forAll(fileNames, i)
    {
        IOdictionary dict
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
        );
        edgeVelocitiesPtrList.append
        (
            new List<timeVelocities>(dict.lookup("edgeVelocities"))
        );
    }
    
    List<timeVelocities> mergedList = 
        mergeEdgeVelocities(edgeVelocitiesPtrList);

    // convert the SLList to lists so that the merged list can be put 
    // into a dictionary
    List<Tuple2<scalar, scalarList> > mergedListConverted(mergedList.size());
    forAll(mergedList, i)
    {
        mergedListConverted[i] = Tuple2<scalar, scalarList>
                                 (
                                    mergedList[i].first(),
                                    scalarList(mergedList[i].second())
                                 );
    }

    IOdictionary outputDict
    (
        IOobject
        (
            outputFile,
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false // do not register in database
        )
    );
    outputDict.set<List<Tuple2<scalar, scalarList> > >("edgeVelocities", mergedListConverted);
    outputDict.regIOobject::write();

    Info<< "Wrote merged edge velocities." << nl
        << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

