/*---------------------------------------------------------------------------*\
 
Application
    testPolygonalTube

Description
    Test the class polygonalTube

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "polygonalTube.H"


int main(int argc, char *argv[])
{

    // Construct edge
    List<point> points(4);
    points[0] = point::zero;
    points[1] = point(1,0,0);
    points[2] = point(1,1,0);
    points[3] = point(1,2,0);

    polygonalEdge edge1(points);

    // Set diameters
    scalarList diameters(3, 0.3);
    diameters[1] = 0.4;
    diameters[2] = 0.5;

    // Construct ellipseAxes
    ellipseAxes startAxes
                (
                    vector(0,0,1),
                    vector(0,1,0),
                    0.15,
                    0.15
                );

    ellipseAxes endAxes
                (
                    vector(1,0,0),
                    vector(0,0,1),
                    0.25,
                    0.25
                );

    // Construct tube
    bool useEffectiveDiameter = false;
    scalar minimalDiameter = 0.4;
    dictionary tubeOptions;
    tubeOptions.add("useEffectiveDiameter", useEffectiveDiameter);
    tubeOptions.add("minimalDiameter", minimalDiameter);

    polygonalTube tube
    (
        edge1, 
        diameters, 
        startAxes, 
        endAxes, 
        tubeOptions
    );

    Info<< "Ostream operator: " << nl << tube << endl;

    // Construct from two points
    polygonalTube tube2(edge1, diameters[0], startAxes, endAxes);
    Info<< "Constructed from one diameter: " << tube2 << endl;

    // Construct without ellipseAxes
    polygonalTube tube3(edge1, diameters, tubeOptions);
    Info<< "Constructed from without ellipseAxes: " << tube3 << endl;

    // Diameter
    Info<< "Diameters: " << tube.segmentDiameters() << endl;
    Info<< "Mean diameter: " << tube.meanDiameter() << endl;
    Info<< "Effective diameter: " << tube.effectiveDiameter() << endl;
    Info<< "Diameter at s = 0.5: " << tube.diameter(0.5) << endl;
    Info<< "Diameter at s = 1.0: " << tube.diameter(1.0) << endl;
    Info<< "Diameter at s = 1.5: " << tube.diameter(1.5) << endl;
    Info<< "Diameter at s = 2.5: " << tube.diameter(2.5) << endl;

    // ellipseAxes
    Info<< "Ellipse axes before s = 0.5:" << nl << tube.ellipseAxesBefore(0.5) << endl;
    Info<< "Ellipse axes after s = 0.5: " << nl << tube.ellipseAxesAfter(0.5) << endl;
    Info<< "Ellipse axes before s = 1.5:" << nl << tube.ellipseAxesBefore(1.5) << endl;
    Info<< "Ellipse axes after s = 1.5: " << nl << tube.ellipseAxesAfter(1.5) << endl;
    Info<< "Ellipse axes before s = 2.5:" << nl << tube.ellipseAxesBefore(2.5) << endl;
    Info<< "Ellipse axes after s = 2.5: " << nl << tube.ellipseAxesAfter(2.5) << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}



