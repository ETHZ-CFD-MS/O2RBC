/*---------------------------------------------------------------------------*\
 
Application
    testRandom

Description
    Test the random number generator from OpenFOAM in parallel.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOdictionary.H"
#include "Random.H"
#include "randomVariable.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

#ifdef USE_RANDOM
        Pout<< "using random " <<endl;
#else
        Pout<< "not using random " <<endl;
#endif

    Random r(1 + Pstream::myProcNo());
    randomVariable rr
    (
        IOdictionary
        (
            IOobject
            (
                "randomDict",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    for(int i=0; i < 10; ++i)
    {
        Pout<< i << ": " << r.scalar01() << endl;
        Pout<< i << ": " << rr.generate() << endl;
    }


    Info<< "End\n" << endl;

    return 0;
}



