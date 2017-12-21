/*---------------------------------------------------------------------------*\
 
Application
    testIOGraph

Description
    Test the class IOgraph.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "IOgraph.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    
    Foam::IOgraph graph
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

    Info<< "nVertices():     " << graph.nVertices() << endl;
    Info<< "nEdges():        " << graph.nEdges() << endl;
    Info<< "vertexIndices(): " << graph.vertexIndices() << endl;
    Info<< "edgeIndices():   " << graph.edgeIndices() << endl;
    Info<< "adjacencyList(): " << graph.adjacencyList() << endl;

    // Info<< "write() :        " << graph.write() << endl;
    Info<< "Ostream operator:" << graph << endl;

    // Test adjacentVertexIndices(...)
    // labelList edgeIndices = graph.edgeIndices();
    // forAll(edgeIndices, i)
    // {
        // labelPair vv = graph.adjacentVertexIndices(edgeIndices[i]);
        // Info << "Edge with index " << edgeIndices[i] << " has neighbor vertices "
             // << "with indices " << vv << endl;
    // }
    // graph.adjacentVertexIndices(5000);

    return 0;
}



