/*---------------------------------------------------------------------------*\
 
Application
    testRBC

Description
    Unit test for class RBC.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "RBC.H"
#include "dissociationCurve.H"
#include "axisymmetricBodyGeometricState.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    dissociationCurve dissCurve(45.0, 2.5);

    axisymmetricBodyGeometricState geomState
    (
        
        vector(0.0, 0.0, 0.0),
        vector(1.0, 0.0, 0.0),
        0.1
    );


    RBC myRBC
    (
        0,
        geomState
    );

    IOdictionary dict
    (
        IOobject
        (
            "RBCDict",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ
        )
    );

    RBC myRBCDict(dict);

    // Info<< "RBC tip: " << myRBC.tip() << endl;
    // Info<< "RBC direction: " << myRBC.direction() << endl;
    Info<< "RBC constructed from components: " << endl;
    Info<< myRBC << endl;
    Info<< "RBC center: " << myRBC.center() << "\n" << endl;


    Info<< "RBC constructed from dictionary: " << endl;
    Info<< myRBCDict << endl;
    Info<< "RBC center: " << myRBCDict.center() << "\n" << endl;

    // TODO: test constructor from Istream
    RBC myRBCIstream(dict.readStream("RBCDict"));

    Info<< "RBC constructed from Istream: " << endl;
    Info<< myRBCIstream << endl;
    Info<< "RBC center: " << myRBCIstream.center() << "\n" << endl;

    return 0;
}


