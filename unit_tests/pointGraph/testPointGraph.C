/*---------------------------------------------------------------------------*\
 
Application
    testPointGraph

Description
    Test the class pointGraph

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "pointGraph.H"


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    
    pointGraph graph
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

    Info<< "\nTest vertexPosition, adjacentVertexIndices and adjacentEdgeIndices:" 
        << nl << endl;

    forAll(graph.externalVertexIndices(), vI)
    {
        label adjacentEI = graph.adjacentEdgeIndices(vI)[0];

        Info<< "Vertex with internal index " << vI
            << " and external index " << graph.externalVertexIndices()[vI]
            << " has convertex internal index " 
            << graph.externalToInternalVertexIndex(graph.externalVertexIndices()[vI])
            << "," << nl << "position " << graph.vertexPosition(vI) << nl
            << " adjacent vertex indices: " << graph.adjacentVertexIndices(vI)
            << " adjacent edge indices:   " << graph.adjacentEdgeIndices(vI) << nl
            << " other adjacent edge index: " << graph.otherAdjacentEdgeIndex(vI, adjacentEI)
            << endl;
    }

    Info<< "\nTest edgeVertexIndices and edgeVertexPositions:"
        << nl << endl;

    forAll(graph.externalEdgeIndices(), eI)
    {
        labelPair   vpi = graph.edgeVertexIndices(eI);
        Pair<point> vpp = graph.edgeVertexPositions(eI);
        Info<< "Edge with internal index " << eI
            << " and external index " << graph.externalEdgeIndices()[eI]
            << " has convertex internal index " 
            << graph.externalToInternalEdgeIndex(graph.externalEdgeIndices()[eI])
            << "," << nl << "neighbor vertices with indices " 
            << vpi.first() << ", " << vpi.second() << endl
            << "and positions " << vpp.first() << ", " << vpp.second() << endl;

    }

    Info<< "\nEnd\n" << endl;

    return 0;
}



