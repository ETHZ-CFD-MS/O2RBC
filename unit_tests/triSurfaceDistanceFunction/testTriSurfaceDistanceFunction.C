/*---------------------------------------------------------------------------*\
 
Application
    testTriSurfaceDistanceFunction

Description
    Unit test for testTriSurfaceDistanceFunction.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "triSurfaceDistanceFunction.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field P_50_ext\n" << endl;
    volScalarField P_50_ext
    (
        IOobject
        (
            "P_50_ext",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    List<scalar> params(3);
    params[0] = 29.3;
    params[1] = 3.0;
    params[2] = 2e-7;

    Info<< "Construct from components..." << endl;
    {
        const scalar threshold = 1e-6;
        const fileName surfName = "RBC_1mm_s.stl";

        triSurfaceDistanceFunction surfDistFcn
        (
            mesh,
            surfName,
            threshold,
            params
        );
    }

    Info<< "Construct from dictionary..." << endl;
    const word dictName("triSurfaceDistanceFunctionDict");

    IOdictionary triSurfaceDistanceFunctionDict
    (
        IOobject
        (
            dictName,
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    triSurfaceDistanceFunction surfDistFcn
    (
        mesh,
        triSurfaceDistanceFunctionDict
    );

    Info<< "Evaluate function..." << endl;
    surfDistFcn.evaluate(P_50_ext);

    runTime += 1.0;

    P_50_ext.write();




    return 0;
}


