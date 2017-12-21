/*---------------------------------------------------------------------------*\
 
Application
    testCylinderGraph

Description
    Test the class cylinderGraph.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "cylinderGraph.H"


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    
    cylinderGraph graph
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

    Info<< "\nTest edgeVertexIndices, edgeVertexPositions, edgeVector"
        << "and edgeLength:" << endl;
    labelList edgeIndices = graph.edgeIndices();
    forAll(edgeIndices, i)
    {
        label eI = edgeIndices[i];
        labelPair   vpi = graph.edgeVertexIndices(eI);
        Pair<point> vpp = graph.edgeVertexPositions(eI);
        vector edgeV1 = graph.edgeVector(eI, vpi.first());
        vector edgeV2 = graph.edgeVector(eI, vpi.second());
        Info<< "Edge with index " << edgeIndices[i] << " has neighbor vertices "
            << "with indices " << vpi.first() << ", " << vpi.second() << endl
            << "and positions " << vpp.first() << ", " << vpp.second() << endl
            << "Edge vector starting from vertex " 
            << vpi.first() << ": " << edgeV1 << endl
            << "Edge vector starting from vertex " 
            << vpi.second() << ": " << edgeV2 << endl
            << "Length = " << graph.edgeLength(eI) << endl;
    }

    Info<< "\nTest vertexDiameter, vertexPosition, adjacentEdges "
        << "and otherAdjacentEdge:" << endl;
    labelList vertexIndices = graph.vertexIndices();
    forAll(vertexIndices, i)
    {
        label vI = vertexIndices[i];
        labelList adjE = graph.adjacentEdges(vI);
        Info<< "Diameter of vertex " << vI << ": " << graph.vertexDiameter(vI) << endl;
        Info<< "Position of vertex " << vI << ": " << graph.vertexPosition(vI) << endl;
        Info<< "Adjacent edges: " << adjE << endl;
        Info<< "Adjacent edge other than " << adjE[0] << ": " 
            << graph.otherAdjacentEdge(vI, adjE[0]) << endl; 
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}



