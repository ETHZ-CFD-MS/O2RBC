/*---------------------------------------------------------------------------*\
 
Application
    testPrimitiveGraph

Description
    Test the class primitiveGraph.
    This code is based on the quick tour found in
    http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/quick_tour.html

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "primitiveGraph.H"

#include <iostream>

struct VProp
{
    int index;
};

struct EProp
{
    int index;
};

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    
    primitiveGraph<VProp, EProp> graph
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

    const labelList& vertexIndices = graph.vertexIndices();
    const labelList& edgeIndices   = graph.edgeIndices();

    Info<< "\nTest edgeVertices(label edgeIndex):" << endl;
    forAll(edgeIndices, i)
    {
        primitiveGraph<VProp, EProp>::VertexD v;
        Pair<primitiveGraph<VProp, EProp>::VertexD> vv 
                                = graph.edgeVertices(edgeIndices[i]);
        v = vv.first();
        std::cout << "Edge with index " << edgeIndices[i] 
                  << " has neighbor vertices with indices " 
                  << int(vv.first()) << ", " << int(vv.second()) 
                  << std::endl;
    }

    Info<< "\nTest adjacentVertices(label vertexIndex):" << endl;
    forAll(vertexIndices, i)
    {
        List<primitiveGraph<VProp, EProp>::VertexD> adjacentList 
            = graph.adjacentVertices(vertexIndices[i]);

        std::cout << "Adjacent vertices to " << vertexIndices[i] << ": ";
        forAll(adjacentList, j)
        {
            std::cout << adjacentList[j] << " ";
        }
        std::cout << std::endl;
    }

    // Test degree
    Info<< "\nTest degree(label vertexIndex):" << endl;
    forAll(vertexIndices, i)
    {
        std::cout << "Degree of vertex " << vertexIndices[i] << ": "
                  << graph.degree(vertexIndices[i]) << std::endl;
    }
    
    Info<< "\nEnd\n" << endl;

    return 0;
}



