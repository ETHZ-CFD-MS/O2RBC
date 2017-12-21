/*---------------------------------------------------------------------------*\
 
Application
    testGeometricOps

Description
    Test the class geometricOps

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "geometricOps.H"


int main(int argc, char *argv[])
{
    vector v1(1,0,0);
    vector v2(0,1,0);
    vector v3(0,0,1);
    vector v4(1,1,0);

    Info<< "normalVector(v1): " << normalVector(v1) << endl;
    Info<< "normalVector(v2): " << normalVector(v2) << endl;
    Info<< "normalVector(v3): " << normalVector(v3) << endl;

    Info<< "normalVector(v1,v2): " << normalVector(v1,v2) << endl;
    Info<< "normalVector(v1,v3): " << normalVector(v1,v3) << endl;
    Info<< "normalVector(v2,v3): " << normalVector(v2,v3) << endl;

    Info<< "bisectionVector(v1,v2): " << bisectionVector(v1,v2) << endl;
    Info<< "bisectionVector(v1,v3): " << bisectionVector(v1,v3) << endl;
    Info<< "bisectionVector(v1,v4): " << bisectionVector(v1,v4) << endl;
    Info<< "bisectionVector(v1,-v1): " << bisectionVector(v1,-v1) << endl;
    Info<< "bisectionVector(v1,v1): " << bisectionVector(v1,v1) << endl;

    Info<< "majorBisectionVector(v1,v2): " << majorBisectionVector(v1,v2) << endl;
    Info<< "majorBisectionVector(v1,v3): " << majorBisectionVector(v1,v3) << endl;
    Info<< "majorBisectionVector(v1,v4): " << majorBisectionVector(v1,v4) << endl;
    Info<< "majorBisectionVector(v1,-v1): " << majorBisectionVector(v1,-v1) << endl;
    Info<< "majorBisectionVector(v1,v1): " << majorBisectionVector(v1,v1) << endl;

    Info<< "bisectionAngle(v1,v2): " << bisectionAngle(v1,v2) << endl;
    Info<< "bisectionAngle(v1,v3): " << bisectionAngle(v1,v3) << endl;
    Info<< "bisectionAngle(v1,v4): " << bisectionAngle(v1,v4) << endl;
    Info<< "bisectionAngle(v1,-v1): " << bisectionAngle(v1,-v1) << endl;
    Info<< "bisectionAngle(v1,v1): " << bisectionAngle(v1,v1) << endl;

    Info<< "majorBisectionAngle(v1,v2): " << majorBisectionAngle(v1,v2) << endl;
    Info<< "majorBisectionAngle(v1,v3): " << majorBisectionAngle(v1,v3) << endl;
    Info<< "majorBisectionAngle(v1,v4): " << majorBisectionAngle(v1,v4) << endl;
    Info<< "majorBisectionAngle(v1,-v1): " << majorBisectionAngle(v1,-v1) << endl;
    Info<< "majorBisectionAngle(v1,v1): " << majorBisectionAngle(v1,v1) << endl;

    Info<< "rotationTensor((v1,v2),(v1,v3)): " 
        << rotationTensor(Pair<vector>(v1,v2), Pair<vector>(v1,v3)) << endl;
    Info<< "rotationTensor((v1,v2),(v2,v3)): " 
        << rotationTensor(Pair<vector>(v1,v2), Pair<vector>(v2,v3)) << endl;

    vector w1(1,1,0);
    Info<< "shearTensor(v1,w1): " << shearTensor(v1, w1) << endl;
    Info<< "shearTensor(v1,-w1): " << shearTensor(v1, -w1) << endl;


    Info<< "\nEnd\n" << endl;

    return 0;
}



