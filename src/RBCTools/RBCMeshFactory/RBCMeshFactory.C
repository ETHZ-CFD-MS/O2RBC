/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "RBCMeshFactory.H"

#include "pointZone.H"
#include "cellZone.H"
#include "faceZone.H"
#include "ZoneMesh.H"

#include "IFstream.H"
#include "deformableBodyGeometricState.H"


// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

namespace Foam
{
    const word RBCMeshFactory::RBCPrefix = "RBC";
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

autoPtr<polyMesh>
Foam::RBCMeshFactory::createRBCMesh
(
    const word& name,
    const Time& runTime,
    const label index
)
{
    // use a read-constructor with
    //    READ_IF_PRESENT: read points, faces, neighbour and owner, but no zones
    //    no parallel communications: this mesh should be created only on one
    //                                processor
    autoPtr<polyMesh> pMesh
    (
        new polyMesh
        (   
            IOobject
            (
                name,
                runTime.caseConstant(),
                runTime,
                IOobject::READ_IF_PRESENT, // does not read zones
                IOobject::NO_WRITE,
                false // do not register in database
            ),
            xferCopy(pointField()),
            xferCopy(faceList()),
            xferCopy(labelList()),
            xferCopy(labelList()),
            false // no parallel comms
        )
    );

    // Add zones. This assumes that the mesh pMesh() does not have any zone. It
    // is up to the user to ensure this (currently, there is no way to remove
    // or rename zones at run time).
    pointZoneMesh& pzm = pMesh->pointZones();
    List<pointZone*> pz;
    pz.append
    (
        new pointZone
        (
            pointZoneName(index),
            identity(pMesh->nPoints()),
            0,
            pzm
        )
    );

    cellZoneMesh& czm = pMesh->cellZones();
    List<cellZone*> cz;
    cz.append
    (
        new cellZone
        (
            cellZoneName(index),
            identity(pMesh->nCells()),
            0,
            czm
        )
    );

    faceZoneMesh& fzm = pMesh->faceZones();
    List<faceZone*> fz;
    fz.append
    (
        new faceZone
        (
            faceZoneName(index),
            identity(pMesh->nFaces()),
            boolList(pMesh->nFaces(), true),
            0,
            fzm
        )
    );

    pMesh->addZones(pz, fz, cz);

    return pMesh;
}

word 
Foam::RBCMeshFactory::pointZoneName(label i)
{
    std::stringstream sstm;
    sstm << RBCPrefix << i << "Points";
    return sstm.str();
}


word 
Foam::RBCMeshFactory::cellZoneName(label i)
{
    std::stringstream sstm;
    sstm << RBCPrefix << i << "Cells";
    return sstm.str();
}


word 
Foam::RBCMeshFactory::faceZoneName(label i)
{
    std::stringstream sstm;
    sstm << RBCPrefix << i << "Faces";
    return sstm.str();
}


autoPtr<deformableBodyGeometricState> 
Foam::RBCMeshFactory::geometricState
(
    const polyMesh& mesh
)
{
    // do not use an IOdictionary to read the dictionary to avoid parallel
    // reading
    IFstream is
    (
        IOobject
        (
            "RBCBaseGeometricState",
            mesh.time().caseConstant(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false // do not register in database
        ).filePath()
    );
    dictionary dict(is);

    return deformableBodyGeometricState::New(dict, mesh);
}


// ************************************************************************* //
