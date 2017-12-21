/*---------------------------------------------------------------------------*\
 
Application
    testEllipseAxes

Description
    Test the class ellipseAxes

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "ellipseAxes.H"


int main(int argc, char *argv[])
{
    ellipseAxes axes1
    (
        vector(1,0,0),
        vector(0,1,0),
        0.3,
        0.5
    );

    ellipseAxes axes2
    (
        // vector(1,0,0), // test orthogonality check
        // vector(0,2,0), // test norm check
        vector(0,1,0), // valid version 
        vector(-1,0,0),
        0.4,
        0.4
    );

    ellipseAxes axes3
    (
        vector(-1,0,0),
        vector(0,-1,0),
        0.3,
        0.5
    );

    ellipseAxes axes4
    (
        vector(Foam::sqrt(2.)/2,Foam::sqrt(2.)/2,0),
        vector(-Foam::sqrt(2.)/2,Foam::sqrt(2.)/2,0),
        0.3,
        0.5
    );

    ellipseAxes axes5
    (
        vector(0,0,1),
        vector(-1,0,0),
        0.3,
        0.5
    );

    Info<< "Ostream operator: axis1 : " << nl << axes1 << endl;
    Info<< "Ostream operator: axis2 : " << nl << axes2 << endl;

    Info<< "Interpolated axes between axes1 and axes2: " 
        << nl << axes1.interpolate(0.5, axes2) << endl;
    Info<< "Interpolated axes between axes1 and axes3: " 
        << nl << axes1.interpolate(0.5, axes3) << endl;
    Info<< "Interpolated axes between axes1 and axes4: " 
        << nl << axes1.interpolate(0.00, axes4)
        << nl << axes1.interpolate(0.25, axes4)
        << nl << axes1.interpolate(0.50, axes4)
        << nl << axes1.interpolate(0.75, axes4)
        << nl << axes1.interpolate(1.00, axes4) << endl;

    Info<< "Interpolated axes between axes1 and axes5: " 
        << nl << axes1.interpolate(0.50, axes5) << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}



