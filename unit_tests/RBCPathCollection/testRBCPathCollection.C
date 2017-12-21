/*---------------------------------------------------------------------------*\
 
Application
    testRBCPathCollection

Description
    Test the class RBCPathCollection

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "RBCPathCollection.H"


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    RBCPathCollection RBCPaths
    (
        IOobject
        (
            "RBCPaths",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    Info<< "\nOstream operator: " << nl << RBCPaths << nl << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info<< "Active paths       : " 
            << RBCPaths.activeIndices() << endl;
        Info<< "Inactive paths     : " 
            << RBCPaths.inactiveIndices() << endl;
        Info<< "New active paths   : " 
            << RBCPaths.newActiveIndices() << endl;
        Info<< "New inactive paths : " 
            << RBCPaths.newInactiveIndices() << nl << endl;
    }

    Info<< "End\n" << endl;
}
    

