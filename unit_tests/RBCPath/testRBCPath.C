/*---------------------------------------------------------------------------*\
 
Application
    testRBCPath

Description
    Test the class RBCPath

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "RBCPath.H"


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    IOdictionary dict
    (
        IOobject
        (
            "sampleRBCPath",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    RBCPath samplePath(dict);

    scalar deltaT = 0.5;

    Info<< "Ostream operator: " << nl << samplePath << nl << endl;

    for(scalar time = -1; time < 4; time += 1)
    {
        scalar oldTime = time - deltaT;

        Info<< "Path active at time " << time << "         : " 
            << samplePath.active(time) << endl;
        Info<< "Path newly active at time " << time << "   : " 
            << samplePath.newActive(time, oldTime) << endl;
        Info<< "Path newly inactive at time " << time << " : " 
            << samplePath.newInactive(time, oldTime) << nl << endl;
    }

}
    

