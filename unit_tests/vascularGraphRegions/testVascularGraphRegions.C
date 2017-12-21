/*---------------------------------------------------------------------------*\
 
Application
    testVascularGraphRegions

Description
    Test the class vascularGraphRegions

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "vascularGraphRegions.H"

#include "circularTubeGraph.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    circularTubeGraph graph
    (
        Foam::IOobject
        (
            "graphDict",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info << graph << endl << endl;

    vascularGraphRegions vesselRegions(mesh, graph, "cylinder");

    vesselRegions.writeMeshes();
    vesselRegions.writeFields();

    Pout<< "\nEnd" << endl;

    return 0;
}



