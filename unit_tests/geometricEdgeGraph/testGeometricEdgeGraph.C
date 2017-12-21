/*---------------------------------------------------------------------------*\
 
Application
    testgeometricEdgeGraph

Description
    Test the class geometricEdgeGraph

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "geometricEdgeGraph.H"


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

    Info<< "Ostream operator: " << graph << endl;

    Info<< "\nTest edgeTangentVectorFromVertex, normalToVertexPlane, "
        << "bisectionVector and bisectionAngle:" 
        << nl << endl;

    forAll(graph.externalVertexIndices(), vI)
    {
        labelList adjacentEdges = graph.adjacentEdgeIndices(vI);

        Info<< "Vertex with internal index " << vI
            << " and external index " << graph.externalVertexIndices()[vI]
            << ": " << endl;

        forAll(adjacentEdges, i)
        {
            label eI = adjacentEdges[i];
            Info<< "Looking at edge " << eI << "..." << endl;

            Info<< "edgeTangentVectorFromVertex: " 
                << graph.edgeTangentVectorFromVertex(vI, eI) << nl
                << "normalToVertexPlane: "
                << graph.normalToVertexPlane(vI, eI) << nl
                << "bisectionVector: "
                << graph.bisectionVector(vI, eI) << nl
                << "bisectionAngle: "
                << graph.bisectionAngle(vI, eI) << nl
                << "majorBisectionVector: "
                << graph.majorBisectionVector(vI, eI) << nl
                << "majorBisectionAngle: "
                << graph.majorBisectionAngle(vI, eI) << nl
                << endl;
        }
    }

    Info<< "\nTest edge accessing:" << nl << endl;

    forAll(graph.externalEdgeIndices(), eI)
    {
        Info<< "Edge with internal index " << eI
            << " and external index " << graph.externalEdgeIndices()[eI]
            << " is: " << nl << graph.edge(eI) << "." << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}



