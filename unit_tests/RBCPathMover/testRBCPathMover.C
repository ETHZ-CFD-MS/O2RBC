/*---------------------------------------------------------------------------*\
 
Application
    testRBCPathMover

Description
    Test the class RBCPathMover and related functionalities:
        - motion of RBCs
        - attaching/detaching RBC paths to/from RBCs
        - interpolation of RBC paths

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "circularTubeGraph.H"
#include "RBC.H"
#include "RBCCollection.H"
#include "RBCPositionTester.H"
#include "RBCPathMover.H"
#include "dissociationCurveHill.H"
#include "regionProperties.H"

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

    dissociationCurveHill dissCurve = dissociationCurveHill(45., 2.5);
    RBCCollection RBCs(mesh, dissCurve);
    RBCPositionTester positionTester(mesh, RBCs);
    RBCPathMover mover
    (
        mesh,
        RBCs,
        positionTester,
        graph
    );

    mover.setInitialPositions();

    Info<< "Finished setting initial positions..." << endl;
    Info<< RBCs << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        mover.moveAll();

        Info<< RBCs << endl;
        Info<< "RBC total mesh volume: " << sum(RBCs.RBCMesh().V()) << nl << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}



