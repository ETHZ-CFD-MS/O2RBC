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
    computeTissueGraphRegions

Description
    Compute the tissue regions closest to each edge of a tube graph.

\*---------------------------------------------------------------------------*/

#include <iomanip>
#include <fstream>

#include "fvCFD.H"

#include "circularTubeGraph.H"
#include "vascularGraphRegions.H"
#include "tissueGraphRegions.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    argList::addNote
    (
        "compute the tissue regions closest to each edge of a tube graph"
    );

    argList::addOption
    (
        "graphDict",
        "word",
        "specify the graph dictionary"
    );

    argList::addOption
    (
        "volumeOutputFile",
        "word",
        "specify the file name for writing tissue volumes"
    );

    word graphDictName = "graphDict";
    if (args.optionFound("graphDict"))
    {
        graphDictName = args["graphDict"];
        Info<< "Graph dictionary: " << graphDictName << endl;
    }

    word volumeOutputFile = "topologicalTissueVolumes.txt";
    if (args.optionFound("volumeOutputFile"))
    {
        volumeOutputFile = args["volumeOutputFile"];
        Info<< "Output file name for tissue volumes: " << volumeOutputFile << endl;
    }

    #include "createTime.H"
    #include "createMesh.H"

    // Create graph
    circularTubeGraph graph
    (
        Foam::IOobject
        (
            graphDictName,
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read PO2 for averaging
    volScalarField PO2
    (
        IOobject
        (
            "PO2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Read the vascular graph regions
    vascularGraphRegions vesselRegions(mesh, graph);

    // Create the tissue graph regions
    tissueGraphRegions tissueRegions(mesh, graph);

    scalarField tissueCellVolumes
    (
        vesselRegions.inTissue().internalField()*mesh.V().field()
    );

    Map<scalar> tissueVolumes = tissueRegions.sumOnTissueRegions(tissueCellVolumes);

    Info<< "Summed tissue region volumes: " << tissueVolumes << nl << endl;

    Info<< "Averaged PO2 on tissue regions: "
        << tissueRegions.averageOnTissueRegions(PO2, tissueCellVolumes) << nl << endl;

    Info<< "Total tissue volume: "
        << gSum(vesselRegions.inTissue().internalField()*mesh.V().field()) << endl;

    // Write the output
    Info<< "Writing fields\n" << endl;
    tissueRegions.writeFields();

    OFstream filestream(volumeOutputFile);
    for 
    (
        Map<scalar>::const_iterator iter = tissueVolumes.begin(); 
        iter != tissueVolumes.cend();
        ++iter
    )
    {
        filestream<< iter.key() << ' ' << *iter << endl;
    }
    Info<< "Wrote tissue volumes to " << volumeOutputFile << "\n" << endl;

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}

