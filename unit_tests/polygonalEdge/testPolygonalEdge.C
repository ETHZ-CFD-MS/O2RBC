/*---------------------------------------------------------------------------*\
 
Application
    testPolygonalEdge

Description
    Test the class polygonalEdge

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "polygonalEdge.H"


int main(int argc, char *argv[])
{

    // Construct from data
    List<point> points(3);
    points[0] = point::zero;
    points[1] = point(1,0,0);
    points[2] = point(1,1,0);

    polygonalEdge edge1(points);

    Info<< "Ostream operator: " << nl << edge1 << endl;

    // Construct from two points
    polygonalEdge edge11(point(0,0,0), point(1,1,1));
    Info<< "Constructed from two points: " << edge11 << endl;

    // Copy constructor
    polygonalEdge edge2(edge1);
    Info<< "Copy of edge1: " << edge2 << endl;

    // operator=
    List<point> points3 = points;
    points3.append(point(2,1,0));
    polygonalEdge edge3(points3);
    edge3 = edge1;
    Info<< "Result of operator= : " << edge3 << endl;

    // length
    Info<< "Length of edge1 = " << edge1.length() << endl;

    // pointPosition
    Info<< "Point position for s = 0.5: " << edge1.pointPosition(0.5) << endl;
    Info<< "Point position for s = 1:   " << edge1.pointPosition(1.0) << endl;
    Info<< "Point position for s = 1.5: " << edge1.pointPosition(1.5) << endl;

    // tangent vector
    Info<< "Tangent vector for s = 0.5: " << edge1.tangentVector(0.5) << endl;
    Info<< "Tangent vector for s = 1:   " << edge1.tangentVector(1.0) << endl;
    Info<< "Tangent vector for s = 1.5: " << edge1.tangentVector(1.5) << endl;

    // segmentEndPoints
    Info<< "Segment end points for s = 0.5: " << edge1.segmentEndPoints(0.5) << endl;
    Info<< "Segment end points for s = 1:   " << edge1.segmentEndPoints(1.0) << endl;
    Info<< "Segment end points for s = 1.5: " << edge1.segmentEndPoints(1.5) << endl;

    // segmentIndex
    Info<< "Segment index for s = 0.5: " << edge1.segmentIndex(0.5) << endl;
    Info<< "Segment index for s = 1:   " << edge1.segmentIndex(1.0) << endl;
    Info<< "Segment index for s = 1.5: " << edge1.segmentIndex(1.5) << endl;

    // segmentLenghts
    Info<< "Segment lengths: " << edge1.segmentLengths() << endl;

    // pathVertexCoords
    Info<< "Path vertex coordinates: " << edge1.pathVertexCoords() << endl;

    // tangentVectorFromPathVertex
    Info<< "Tangent vector from path vertex 0 in positive direction: " 
        << edge1.tangentVectorFromPathVertex(0, true) << endl;
    Info<< "Tangent vector from path vertex 1 in negative direction: " 
        << edge1.tangentVectorFromPathVertex(1, false) << endl;
    Info<< "Tangent vector from path vertex 1 in positive direction: " 
        << edge1.tangentVectorFromPathVertex(1, true) << endl;
    Info<< "Tangent vector from path vertex 2 in negative direction: " 
        << edge1.tangentVectorFromPathVertex(2, false) << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}



