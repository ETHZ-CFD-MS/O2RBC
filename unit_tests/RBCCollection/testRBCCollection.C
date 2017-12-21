/*---------------------------------------------------------------------------*\
 
Application
    testRBC

Description
    Unit test for class RBC.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "RBC.H"
#include "RBCCollection.H"
#include "dissociationCurveHill.H"
#include "axisymmetricBodyGeometricState.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    dissociationCurveHill dissCurve = dissociationCurveHill(45., 2.5);
    Info<< "Bounding box of the Eulerian mesh: " << mesh.bounds() << endl;

    // construct RBCCollection from empty RBCCollection file
    RBCCollection RBCs
    (
        mesh,
        dissCurve
    );

    // TEST: RBC insertion
    runTime++;

    RBCs.insert(1, axisymmetricBodyGeometricState(point(   -2,1,1), vector(1, 0,0), 3e-6, 20e-6), 90);
    RBCs.insert(3, axisymmetricBodyGeometricState(point(-1e-5,0,0), vector(1, 0,0), 3e-6, 20e-6), 90);
    RBCs.insert(5, axisymmetricBodyGeometricState(point( 2e-5,0,0), vector(1,-1,0), 3e-6, 20e-6), 90);
    RBCs.insert(7, axisymmetricBodyGeometricState(point( 4e-5,0,0), vector(1, 0,0), 3e-6, 20e-6), 90);
    RBCs.updateMeshes();
    runTime.write();

    // TEST: RBC deletion
    runTime++;
    RBCs.remove(5);
    RBCs.updateMeshes();
    runTime.write();

    // TEST: RBC insertion
    runTime++;
    RBCs.insert(9, axisymmetricBodyGeometricState(point(6e-5,0,0), vector(0,0,1), 3e-6, 20e-6), 90);
    RBCs.updateMeshes();
    runTime.write();

    // TEST: RBC deletion
    runTime++;
    RBCs.remove(7);
    RBCs.updateMeshes();
    runTime.write();

    // TEST: simultaneous RBC insertion and deletion
    runTime++;
    RBCs.remove(9);
    RBCs.insert(11, axisymmetricBodyGeometricState(point(3e-5,0,0), vector(1,0,0), 3e-6, 20e-6), 90);
    RBCs.updateMeshes();
    runTime.write();

    // TEST: RBC motion
    runTime++;
    RBCs.prepareMotion(1, axisymmetricBodyGeometricState(point(-1,0,0), vector(1,0,0), 3e-6, 20e-6));
    RBCs.prepareMotion(3, axisymmetricBodyGeometricState(point(2e-5,0,0), vector(1,0,0), 3e-6, 20e-6));
    RBCs.prepareMotion(11, axisymmetricBodyGeometricState(point(5e-5,0,0), vector(1,0,0), 3e-6, 20e-6));
    RBCs.applyMotion();
    runTime.write();

    // TEST: RBC motion (RBC1 enters the domain!)
    runTime++;
    RBCs.prepareMotion(1, axisymmetricBodyGeometricState(point(0,0,0), vector(1,0,0), 3e-6, 20e-6));
    RBCs.prepareMotion(3, axisymmetricBodyGeometricState(point(3e-5,1e-5,0), vector(1,0,0), 3e-6, 20e-6));
    RBCs.prepareMotion(11, axisymmetricBodyGeometricState(point(6e-5,2e-5,0), vector(1,0,0), 3e-6, 20e-6));
    RBCs.applyMotion();
    runTime.write();

    // TEST: RBC motion (RBC11 leaves the domain!)
    runTime++;
    RBCs.prepareMotion(1, axisymmetricBodyGeometricState(point(1e-5,0,0), vector(1,0.1,0), 3e-6, 20e-6));
    RBCs.prepareMotion(3, axisymmetricBodyGeometricState(point(4e-5,1e-5,0), vector(1,0.1,0), 3e-6, 20e-6));
    RBCs.prepareMotion(11, axisymmetricBodyGeometricState(point(2,2e-5,0), vector(1,0,0), 3e-6, 20e-6));
    RBCs.applyMotion();
    runTime.write();

    // TEST: construct object from nonempty RBCCollection file (use previous
    // time step)
    RBCCollection RBCs2
    (
        mesh,
        dissCurve
    );
    runTime++;
    Info<< RBCs2;
    runTime.write();



    return 0;
}


