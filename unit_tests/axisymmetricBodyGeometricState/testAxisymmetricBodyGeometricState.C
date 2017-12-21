/*---------------------------------------------------------------------------*\
 
Application
    testAxisymmetricBodyGeometricState

Description
    Test the class axisymmetricBodyGeometricState

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "axisymmetricBodyGeometricState.H"
#include "IOdictionary.H"


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    
    IOdictionary baseDict
    (
        Foam::IOobject
        (
            "baseGeometricState",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    IOdictionary dictForMeshConstruction
    (
        Foam::IOobject
        (
            "baseGeometricState",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    polyMesh referenceMesh
    (
        Foam::IOobject
        (
            "referenceMesh",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    IOdictionary newDict
    (
        Foam::IOobject
        (
            "newGeometricState",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    axisymmetricBodyGeometricState baseGeomState(baseDict);
    axisymmetricBodyGeometricState baseGeomStateFromMesh(baseDict, referenceMesh);
    axisymmetricBodyGeometricState newGeomState(newDict);

    Info<< "Bounding box of base state: " << baseGeomState.bounds() << endl;
    Info<< "Bounding box of mesh state: " << baseGeomStateFromMesh.bounds() << endl;
    Info<< "Bounding box of new  state: " << newGeomState.bounds() << endl;

    Info<< endl;
    
    Info<< "Ostream operator: " << nl 
        << "base state:" << nl << baseGeomState << nl
        << "mesh state:" << nl << baseGeomStateFromMesh << nl
        << "new state: " << nl << newGeomState << endl;

    pointField pf(4);
    pf[0] = vector::zero;
    pf[1] = vector(1.0, 0, 0);
    pf[2] = vector(0.0, 0.1, 0);
    pf[3] = vector(0.0, 0.0, 0.1);

    Info<< "Base points: " << pf << endl;
    Info<< "Transformed points: " << baseGeomState.transformPoints(pf, newGeomState, true) << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}



