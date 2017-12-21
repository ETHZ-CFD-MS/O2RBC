/*---------------------------------------------------------------------------*\
 
Application
    testTubeMapperVisitor

Description
    Test the class tubeMapperVisitor

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "tubeMapperVisitor.H"


int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"

    // Construct edge
    List<point> points(7);
    points[0] = point::zero;
    points[1] = point(1,0,0);
    points[2] = point(1,1,0);
    points[3] = point(1,1,1);
    points[4] = point(0,1,1);
    points[5] = point(0,0,1);
    points[6] = point(0,0,0);

    polygonalEdge edge1(points);

    // Set diameters
    scalarList diameters(6, 0.4);

    // Construct ellipseAxes
    ellipseAxes startAxes
                (
                    vector(0,0,1),
                    vector(0,1,0),
                    0.2,
                    0.2
                );

    ellipseAxes endAxes
                (
                    vector(0,1,0),
                    vector(1,0,0),
                    0.2,
                    0.2
                );

    // Construct tube
    polygonalTube tube(edge1, diameters, startAxes, endAxes);
    Info<< "polygonalTube used for testing: " << nl << tube << endl;

    // Construct tubeMapperVisitor
    tubeMapperVisitor tubeMapper
    (
        IOobject
        (
            "cylinder",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ
        )
    );

    Info<< endl;

    // visit
    tubeMapper.visitPolygonalTube(tube);

    // transformPoint
    Info<<           "transformPoint(point(0,1,0), tube, 0): " 
        << tubeMapper.transformPoint(point(0,1,0), tube, 0) << endl;
    Info<<           "transformPoint(point(1,1,0), tube, 0): " 
        << tubeMapper.transformPoint(point(1,1,0), tube, 0) << endl;
    Info<<           "transformPoint(point(2,1,0), tube, 0): " 
        << tubeMapper.transformPoint(point(2,1,0), tube, 0) << endl;
    Info<<           "transformPoint(point(2,1,0), tube, 1): " 
        << tubeMapper.transformPoint(point(2,1,0), tube, 1) << endl;
    Info<<           "transformPoint(point(3,1,0), tube, 1): " 
        << tubeMapper.transformPoint(point(3,1,0), tube, 1) << endl;
    Info<<           "transformPoint(point(4,1,0), tube, 1): " 
        << tubeMapper.transformPoint(point(4,1,0), tube, 1) << endl;
    Info<<           "transformPoint(point(4,1,0), tube, 2): " 
        << tubeMapper.transformPoint(point(4,1,0), tube, 2) << endl;
    Info<<           "transformPoint(point(5,1,0), tube, 2): " 
        << tubeMapper.transformPoint(point(5,1,0), tube, 2) << endl;
    Info<<           "transformPoint(point(6,1,0), tube, 2): " 
        << tubeMapper.transformPoint(point(6,1,0), tube, 2) << endl;

    Info<< endl;

    // rotationTensorFromSegmentToAxis
    Info<< "rotationTensorFromSegmentToAxis(tube, 0): "
        << tubeMapper.rotationTensorFromSegmentToAxis(tube, 0) << endl;
    Info<< "rotationTensorFromSegmentToAxis(tube, 1): "
        << tubeMapper.rotationTensorFromSegmentToAxis(tube, 1) << endl;
    Info<< "rotationTensorFromSegmentToAxis(tube, 2): "
        << tubeMapper.rotationTensorFromSegmentToAxis(tube, 2) << endl;

    // cylinderPointIDsOnSegment
    // Info<< "cylinderPointIDsOnSegment(tube, 0): "
        // << tubeMapper.cylinderPointIDsOnSegment(tube, 0) << endl;
    // Info<< "cylinderPointIDsOnSegment(tube, 1): "
        // << tubeMapper.cylinderPointIDsOnSegment(tube, 1) << endl;

    autoPtr<polyMesh> pMappedMesh = tubeMapper.mappedMesh(0);
    pMappedMesh->write();

    autoPtr<volScalarField> pAxisCoord = tubeMapper.normalizedCylinderAxisCoord();
    pAxisCoord->write();

    Info<< "\nEnd\n" << endl;

    return 0;
}



