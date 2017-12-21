/*---------------------------------------------------------------------------*\
 
Application
    testStraightEdge

Description
    Test the class straightEdge

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "straightEdge.H"


int main(int argc, char *argv[])
{

    // Construct from data
    straightEdge edge1(point(0,0,0), point(1,1,1));

    Info<< "Ostream operator: " << nl << edge1 << endl;

    // Copy constructor
    straightEdge edge2(edge1);
    Info<< "Copy of edge1 = " << edge2 << endl;

    // operator=
    straightEdge edge3(point::zero, point::zero);
    edge3 = edge1;
    Info<< "Result of operator= : " << edge3 << endl;

    // length
    Info<< "Length of edge1 = " << edge1.length() << endl;

    // pointPosition
    Info<< "Point position for s = 1: " << edge1.pointPosition(1.0) << endl;

    // tangent vector
    Info<< "Tangent vector for s = 1: " << edge1.tangentVector(1.0) << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}



