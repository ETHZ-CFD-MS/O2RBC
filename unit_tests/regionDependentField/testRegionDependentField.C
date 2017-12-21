/*---------------------------------------------------------------------------*\
 
Application
    testCylinderGraph

Description
    Test the class cylinderGraph

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "cylinderGraph.H"
#include "vascularGraphRegions.H"
#include "regionDependentField.H"


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    cylinderGraph graph
    (
        Foam::IOobject
        (
            "cylinderGraph",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info << graph << endl << endl;

    vascularGraphRegions vesselRegions(mesh, graph, "cylinder");

    regionDependentField alpha
    (
        IOobject(
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        vesselRegions
    );

    alpha.write();


    return 0;
}



