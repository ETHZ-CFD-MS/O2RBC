/*---------------------------------------------------------------------------*\
 
Application
    testGraphCoordinateInterpolation

Description
    Test the class graphCoordinateInterpolation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "geometricEdgeGraph.H"
#include "graphCoordinate.H"
#include "graphCoordinateInterpolation.H"


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    
    geometricEdgeGraph graph
    (
        Foam::IOobject
        (
            "graphDict",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    graphCoordinateInterpolation graphInterp(graph);

    graphCoordinate gc1(0, 0.5);
    graphCoordinate gc2(0, 1);
    graphCoordinate gc3(1, 0.5);
    graphCoordinate gc4(1, 1);
    graphCoordinate gc5(2, 0.5);
    graphCoordinate gc6(2, 1);

    Info<< "gc2 between gc1 and gc3: " << graphInterp.isBetweenUpToNeighbourEdge(gc2, gc1, gc3) << endl;
    Info<< "gc2 between gc3 and gc1: " << graphInterp.isBetweenUpToNeighbourEdge(gc2, gc3, gc1) << endl;
    Info<< "gc1 between gc2 and gc3: " << graphInterp.isBetweenUpToNeighbourEdge(gc1, gc2, gc3) << endl;
    Info<< "gc2 between gc4 and gc1: " << graphInterp.isBetweenUpToNeighbourEdge(gc2, gc4, gc1) << endl;
    Info<< "gc2 between gc2 and gc1: " << graphInterp.isBetweenUpToNeighbourEdge(gc2, gc2, gc1) << endl;
    Info<< "gc2 between gc3 and gc4: " << graphInterp.isBetweenUpToNeighbourEdge(gc2, gc3, gc4) << endl;
    Info<< "gc2 between gc3 and gc5: " << graphInterp.isBetweenUpToNeighbourEdge(gc2, gc3, gc5) << endl;
    Info<< "gc5 between gc2 and gc6: " << graphInterp.isBetweenUpToNeighbourEdge(gc5, gc2, gc6) << endl;

    scalarList times(5);
    times[0] = 0; times[1] = 1; times[2] = 2; times[3] = 3; times[4] = 4;
    labelList edges(5);
    edges[0] = 0;     edges[1] = 0;   edges[2] = 1;     edges[3] = 1;   edges[4] = 3;
    scalarList sCoords(5);
    // different possibilities to try out different edge orientations in
    // graphDict
    // sCoords[0] = 0.5; sCoords[1] = 1.0; sCoords[2] = 0.5; sCoords[3] = 1.0; sCoords[4] = 0.5;
    // sCoords[0] = 0.5; sCoords[1] = 1.0; sCoords[2] = 1.0; sCoords[3] = 0.5; sCoords[4] = 0.5;
    sCoords[0] = 1.0; sCoords[1] = 0.5; sCoords[2] = 0.5; sCoords[3] = 1.0; sCoords[4] = 0.5;
    // sCoords[0] = 1.0; sCoords[1] = 0.5; sCoords[2] = 1.0; sCoords[3] = 0.5; sCoords[4] = 0.5;

    Info<< "time at edge 0, sCoord = 0.25: " << 
        graphInterp.computeGraphCoordinatePassingTime(times, edges, sCoords, graphCoordinate(0, 0.25)) << endl;
    // Info<< "time at edge 0, sCoord = 1.49: " << 
        // graphInterp.computeGraphCoordinatePassingTime(times, edges, sCoords, graphCoordinate(0, 1.49)) << endl;
    Info<< "time at edge 1, sCoord = 0.0: " << 
        graphInterp.computeGraphCoordinatePassingTime(times, edges, sCoords, graphCoordinate(1, 0.0)) << endl;
    Info<< "time at edge 1, sCoord = 0.5: " << 
        graphInterp.computeGraphCoordinatePassingTime(times, edges, sCoords, graphCoordinate(1, 0.5)) << endl;
    Info<< "time at edge 1, sCoord = 1.1: " << 
        graphInterp.computeGraphCoordinatePassingTime(times, edges, sCoords, graphCoordinate(1, 1.1)) << endl;
    Info<< "time at edge 1, sCoord = 1.5: " << 
        graphInterp.computeGraphCoordinatePassingTime(times, edges, sCoords, graphCoordinate(1, 1.5)) << endl;
    Info<< "time at edge 3, sCoord = 0.49: " << 
        graphInterp.computeGraphCoordinatePassingTime(times, edges, sCoords, graphCoordinate(3, 0.49)) << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}



