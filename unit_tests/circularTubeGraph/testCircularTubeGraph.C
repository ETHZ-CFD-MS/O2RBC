/*---------------------------------------------------------------------------*\
 
Application
    testCircularTubeGraph

Description
    Test the class circularTubeGraph

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "circularTubeGraph.H"


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    
    circularTubeGraph graph
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

    forAll(graph.externalVertexIndices(), vI)
    {
        labelList adjacentEdges = graph.adjacentEdgeIndices(vI);

        Info<< "Vertex with internal index " << vI
            << " and external index " << graph.externalVertexIndices()[vI]
            << ": " << endl;
        Info<< "vertex diameter: " << graph.vertexDiameter(vI) << nl
            << "vertex radius:   " << graph.vertexRadius(vI) << endl;

        forAll(adjacentEdges, i)
        {
            label eI = adjacentEdges[i];
            Info<< "Looking at edge " << eI << "..." << endl;

            Info<< "edgeDiameterAtVertex: " 
                << graph.edgeDiameterAtVertex(vI, eI) << nl
                << "edgeRadiusAtVertex: "
                << graph.edgeRadiusAtVertex(vI, eI) << nl
                << "vertexTubeAxes: "
                << graph.vertexTubeAxes(vI, eI) << nl
                << endl;
        }
    }

    Info<< "\nTest tube accessing:" << nl << endl;

    forAll(graph.externalEdgeIndices(), eI)
    {
        Info<< "Edge with internal index " << eI
            << " and external index " << graph.externalEdgeIndices()[eI]
            << " has tube: " << nl << graph.tube(eI) << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}



