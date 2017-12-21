// Compute the field xiRBC, which is the area ratio in Eulerian grid cells
// of the area potentially occupied by flowing RBCs, and the grid cell 
// cross-sectional area.
// For those cell at the upper boundary of the RBC, this ratio is
// obtained by dividing the areas of the trapezoidal shapes of the cross
// sectional areas.
// This version works only for cylindrical RBCs.

scalar yMaxRBC = RBCs.RBCMesh().bounds().max().y();

const labelListList& cellPoints = mesh.cellPoints();
const pointField& points = mesh.points();

forAll(mesh.C(), cI)
{
    // compute maximum y-coordinate of cells
    const labelList& pointIndices = cellPoints[cI];
    List<point> cellP(pointIndices.size());
    forAll(cellP, pI)
    {
        cellP[pI] = points[pointIndices[pI]];
    }
    boundBox cellBounds(cellP, /* doReduce */ false);
    scalar yMinCell = cellBounds.min().y();
    scalar yMaxCell = cellBounds.max().y();
    if (yMaxCell < yMaxRBC)
    {
        xiRBC[cI] = 1.0;
    }
    else if (yMinCell > yMaxRBC)
    {
        xiRBC[cI] = 0.0;
    }
    else
    {
        xiRBC[cI] = (yMinCell + yMaxRBC)/(yMinCell + yMaxCell)
                  * (yMaxRBC - yMinCell)/(yMaxCell - yMinCell);
    } 
}

