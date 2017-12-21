Info<< max(phi) << endl;
Info<< "Upwind cells: " << endl;
cellList cells = mesh.cells();
forAll(cells, cI)
{
    cell cCell = cells[cI];
    // find the cell face with the maximum flux
    scalar maxFlux = 0.0;
    label faceIMaxFlux = 0;
    forAll(cCell, j)
    {
        label faceI = cCell[j];
        // TODO: deal with the case when faceI is a boundary face.
        // This line currently causes an out of bound error.
        if (mag(phi[faceI]) > max(1e-20, maxFlux))
        {
            maxFlux = mag(phi[faceI]);
            faceIMaxFlux = faceI;
        }
    }
    if (maxFlux > 1e-20)
    {
        // Info<< "maxFlux = " << maxFlux << ", faceIMaxFlux = " << faceIMaxFlux << endl;
        // set the upwind cell to be the opposite cell
        upwindCells[cI] = (cI == mesh.owner()[faceIMaxFlux] ?
                           mesh.faceNeighbour()[faceIMaxFlux]:
                           mesh.faceOwner()[faceIMaxFlux]);
        // Info<< upwindCells[cI] << " ----> " << cI << endl;
    }
    else
    {
        upwindCells[cI] = cI;
    }
}


