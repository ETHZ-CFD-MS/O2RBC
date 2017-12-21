/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "RBCCollection.H"
#include "RBCMeshFactory.H"

#include "dissociationCurve.H"
#include "deformableBodyGeometricState.H"
#include "geometricOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RBCCollection, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void
RBCCollection::addMesh(const RBC& cRBC)
{
    label RBCIdx = cRBC.index();
    if (debug)
    {
        Pout<< "Foam::RBCCollection::addMesh(const RBC&):"
            << " Adding mesh for RBC " << RBCIdx << endl;
    }

    RBCMesh_.addMesh
    (
        RBCMeshFactory::createRBCMesh
        (
            cRBC.sampleMeshName(),
            RBCMesh_.time(),
            RBCIdx
        )
    );

    RBCMesh_.merge();
}


void
RBCCollection::addMeshIfInDomain(const RBC& cRBC)
{
    bool createMesh = false;
    label RBCIdx = cRBC.index();

    // test whether a mesh should be added
    if (Pstream::parRun())
    {
        const label procI = geometricRBCProc(RBCIdx);
        if (procI == Pstream::myProcNo() && !(RBCHasNonEmptyZone(RBCIdx)))
        {
            createMesh = true;
            if (debug)
            {
                Pout<< "Foam::RBCCollection::addMeshIfInDomain(const label):"
                    << " Going to create mesh for RBC " << RBCIdx << endl;
            }
        }
    }
    else
    {
        const bool inDomain = RBCInDomain(RBCIdx);
        if (inDomain && !(RBCHasNonEmptyZone(RBCIdx)))
        {
            createMesh = true;
        }
    }

    // add mesh zones if required
    if (createMesh)
    {
        addMesh(cRBC);
        addedList_.append(RBCIdx);
    }
}


void
RBCCollection::removeMesh(const RBC& cRBC)
{
    label RBCIdx = cRBC.index();
    if (debug)
    {
        Pout<< "Foam::RBCCollection::removeMesh(const label):"
            << " Removing mesh zones for RBC " << RBCIdx << endl;
    }

    const word czName = RBCMeshFactory:: cellZoneName(RBCIdx);
    const word pzName = RBCMeshFactory::pointZoneName(RBCIdx);
    const word fzName = RBCMeshFactory:: faceZoneName(RBCIdx);

    RBCMesh_.removeMesh(czName, pzName, fzName);
}


void
RBCCollection::moveMeshes()
{
    pointField meshPoints = RBCMesh_.points();

    forAllConstIter(Map<autoPtr<RBC> >, RBCs_, iter)
    {
        const RBC& cRBC = *iter;
        const label RBCIdx = cRBC.index();

        if (RBCHasNonEmptyZone(RBCIdx))
        {
            pointField zonePoints = bluePrintMesh(cRBC.sampleMeshName()).points();
            const deformableBodyGeometricState& geomState = 
                geometricState(cRBC.sampleMeshName());
            zonePoints = geomState.transformPoints
            (
                zonePoints,
                cRBC.geometricState(),
                /* exactDiameter */ false
            );

            // replace corresponding values in meshPoints
            const labelList& pIdxList = RBCMesh_.pointZones()[RBCZoneID(RBCIdx)];

            forAll(pIdxList, i)
            {
                meshPoints[pIdxList[i]] = zonePoints[i];
            }
        }
    }

    // final call to movePoints
    RBCMesh_.fvMesh::movePoints(meshPoints);
}


void
RBCCollection::initializeAddedFields()
{
    // set constant PO2 and Hb on the added zones
    forAll(addedList_, i)
    {
        const label RBCIdx = addedList_[i];
        const scalar PO2Init = (*this)[RBCIdx].inletPO2();
        
        if (debug)
        {
            Pout<< "Foam::RBCCollection::initializeAddedFields():"
                << " Setting RBC PO2 = " << PO2Init 
                << " for RBC " << RBCIdx << endl;
        }

        setRBCVolumeFraction(RBCIdx);
        setConstantPO2(RBCIdx, PO2Init);
        setFieldOldTime(RBCIdx);
        setVolumeOldTime(RBCIdx);
    }
}


void
RBCCollection::setRBCVolumeFraction(const label RBCIdx)
{
    RBCMesh_.setFieldOnZone(in_RBC_, 1.0, RBCZoneID(RBCIdx));
}


void
RBCCollection::setFieldOldTime(const label RBCIdx)
{
    RBCMesh_.setFieldOnZone
    (
        Hb_.oldTime(), 
        RBCMesh_.getFieldOnZone(Hb_, RBCZoneID(RBCIdx)),
        RBCZoneID(RBCIdx)
    );
}


void
RBCCollection::setVolumeOldTime(const label RBCIdx)
{
    if (RBCMesh_.time().timeIndex() > 0)
    {
        volScalarField::Internal&       V0 = RBCMesh_.setV0();
        const volScalarField::Internal& V  = RBCMesh_.V();

        RBCMesh_.setFieldOnZone
        (
            V0,
            RBCMesh_.getFieldOnZone(V, RBCZoneID(RBCIdx)),
            RBCZoneID(RBCIdx)
        );
    }
}


void
RBCCollection::sendFields()
{
    if (Pstream::parRun())
    {
        // Send field data.
        // Each entry of sendList_ contains the RBC indices that are to be sent
        // to the corresponding processor index. A message consists the field data 
        // of the sent RBCs. Only one OPstream is used for each destination 
        // processor, even if multiple RBCs are sent.
        forAll(sendList_, procI)
        {
            const labelList& RBCList = sendList_[procI];
            if (RBCList.size() > 0)
            {
                OPstream toNext
                (
                    Pstream::blocking,
                    procI,
                    /* bufSize */ 0,
                    /* tag */ 2
                );
                forAll(RBCList, i)
                {
                    const label RBCIdx = RBCList[i];

                    if (debug)
                    {
                        Pout<< "Foam::RBCCollection::sendFields(): "
                            << "Sending fields of RBC " << RBCIdx 
                            << " to proc " << procI << endl;
                    }

                    const scalarField f(RBCMesh_.getFieldOnZone(Hb_, RBCZoneID(RBCIdx)));
                    toNext << f;
                }
            }
        }
    }
}


void
RBCCollection::receiveFields()
{
    if (Pstream::parRun())
    {
        // Receive field data.
        // Each entry of recvList_ contain the RBC indices that are to be received
        // by the current processor. Only one message per sending processor is
        // communicated. Each of these message consists of the local field data 
        // on the sending processor.
        forAll(recvList_, procI)
        {
            const labelList& RBCList = recvList_[procI];
            if (RBCList.size() > 0)
            {
                IPstream fromPrevious
                (
                    Pstream::blocking,
                    procI,
                    /* bufSize */ 0,
                    /* tag */ 2
                );
                forAll(RBCList, i)
                {
                    const label RBCIdx = RBCList[i];

                    if (debug)
                    {
                        Pout<< "Foam::RBCCollection::receiveFields(): "
                            << "Receiving fields of RBC " << RBCIdx
                            << " from processor " << procI << endl;
                    }

                    scalarField f;
                    fromPrevious >> f;

                    // insert received field into the field defined on RBCMesh_
                    RBCMesh_.setFieldOnZone(Hb_, f, RBCZoneID(RBCIdx));
                }
            }
        }
    }
}


void
RBCCollection::setOldReceivedFields()
{
    forAll(recvList_, procI)
    {
        labelList RBCList = recvList_[procI];
        forAll(RBCList, i)
        {
            const label RBCIdx = RBCList[i];
            setFieldOldTime(RBCIdx);
            setVolumeOldTime(RBCIdx);
        }
    }
}


bool
RBCCollection::RBCInDomain(const label RBCIdx) const
{
    const RBC& cRBC = (*this)[RBCIdx];
    return deformableBodyInDomain(cRBC.geometricState());
}


bool
RBCCollection::deformableBodyInDomain(const deformableBodyGeometricState& geomState) const
{
    const polyMesh& mesh = procMeshInfo_.mesh();
    if (mesh.bounds().contains(geomState.center())
     || mesh.bounds().overlaps(geomState.bounds()))
    {
        return true;
    }
    return false;
}


label
RBCCollection::geometricRBCProc(const label RBCIdx) const
{
    const RBC& cRBC = (*this)[RBCIdx];
    return deformableBodyProc(cRBC.geometricState());
}


label
RBCCollection::deformableBodyProc(const deformableBodyGeometricState& geomState) const
{
    // 1. Test the center
    label centerProcI = procMeshInfo_.findProcNo(geomState.center());
    if (centerProcI >= 0)
    {
        return centerProcI;
    }

    // 2. Test the bounding box
    labelList overlappingProcIDs = 
                procMeshInfo_.findOverlappingProcIDs(geomState.bounds());
    if (overlappingProcIDs.size() > 0)
    {
        return overlappingProcIDs[0];
    }

    return -1;
}


void
Foam::RBCCollection::addMeshToTable
(
    const word& sampleMeshName,
    autoPtr<polyMesh> meshPtr
)
{
    if (!bluePrintMeshTable_.found(sampleMeshName))
    {
        bluePrintMeshTable_.set(sampleMeshName, meshPtr);
        geomStateTable_.set
        (
            sampleMeshName, 
            RBCMeshFactory::geometricState(meshPtr())
        );
    }
}


const Foam::polyMesh&
Foam::RBCCollection::bluePrintMesh(const word& sampleMeshName) const
{
    if (!bluePrintMeshTable_.found(sampleMeshName))
    {
        autoPtr<polyMesh> meshPtr =
            RBCMeshFactory::createRBCMesh
            (
                sampleMeshName,
                RBCMesh_.time(),
                0
            );
        bluePrintMeshTable_.set
        (
            sampleMeshName, 
            autoPtr<polyMesh>
            (
                new polyMesh
                (   
                    IOobject
                    (
                        sampleMeshName,
                        RBCMesh_.time().caseConstant(),
                        RBCMesh_.time(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false // do not register in database
                    ),
                    xferCopy(meshPtr->points()),
                    xferCopy(meshPtr->faces()),
                    xferCopy(meshPtr->faceOwner()),
                    xferCopy(meshPtr->faceNeighbour()),
                    false // no parallel comms
                )
            )
        );
        geomStateTable_.set
        (
            sampleMeshName, 
            RBCMeshFactory::geometricState(meshPtr())
        );
    }
    return bluePrintMeshTable_[sampleMeshName]();
}


const Foam::deformableBodyGeometricState&
Foam::RBCCollection::geometricState(const word& sampleMeshName) const
{
    return geomStateTable_[sampleMeshName];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBCCollection::RBCCollection
(
    const fvMesh& mesh,
    const dissociationCurve& dissociationCurve
)
:
    regIOobject
    (
        IOobject
        (
            "RBCCollection",
            mesh.time().timeName(),
            "uniform",
            mesh.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    RBCMesh_
    (
        IOobject
        (
            "RBC",
            mesh.time().timeName(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    ),
    Hb_
    (
        IOobject
        (
            "Hb",
            RBCMesh_.time().timeName(),
            RBCMesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        RBCMesh_
    ),
    PO2_
    (
        IOobject
        (
            "PO2_RBC",
            RBCMesh_.time().timeName(),
            RBCMesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        RBCMesh_
    ),
    in_RBC_
    (
        IOobject
        (
            "in_RBC",
            RBCMesh_.time().timeName(),
            RBCMesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        RBCMesh_,
        dimless
    ),
    dissociationCurve_(dissociationCurve),
    RBCs_(),
    procMeshInfo_(mesh),
    sendList_(Pstream::nProcs()),
    recvList_(Pstream::nProcs()),
    addedList_(),
    bluePrintMeshTable_(),
    geomStateTable_()
{
    in_RBC_ = 1.0;

     // read RBC objects from the RBCCollection dictionary
     const IOdictionary dict
     (
         IOobject
         (
             this->name(),
             this->time().timeName(),
             this->db(),
             IOobject::NO_READ,
             IOobject::NO_WRITE,
             false
         ),
         this->readStream(typeName)
     );
 
     this->close();
 
    // create RBCs
    wordList dictKeys = dict.toc();
    forAll(dictKeys, i)
    {
        const word& key = dictKeys[i];
        if (dict.isDict(key))
        {
            const dictionary& RBCDict = dict.subDict(key); 
            const label RBCIdx = readLabel(RBCDict.lookup("index"));
            RBCs_.insert
            (
                RBCIdx,
                autoPtr<RBC>
                (
                    new RBC(RBCDict)
                )
            );
        }
    }

    // create mesh zones for those RBCs that are inside the domain
    forAllConstIter(Map<autoPtr<RBC> >, RBCs_, iter)
    {
        addMeshIfInDomain(*iter);
    }

    updateMeshes();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBCCollection::~RBCCollection()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

volScalarField&
RBCCollection::getField(
    const word& fieldName
)
{
    if (fieldName == "PO2")
    {
        return PO2_;
    }
    else if (fieldName == "Hb")
    {
        return Hb_;
    }
    else
    {
        FatalErrorIn
        (
            "RBC::getField(const word& fieldName)"
        )   << "  Non constant access is only given for ""PO2"" or " 
            << """Hb""." << nl
            << abort(FatalError);
    }
    // just here to avoid warning
    return PO2_;
}


label
RBCCollection::RBCZoneID(const label RBCIdx) const
{
    if (RBCHasNonEmptyZone(RBCIdx))
    {
        const word czName = RBCMeshFactory::cellZoneName(RBCIdx);
        return RBCMesh_.cellZones().findZoneID(czName);
    }
    else
    {
        FatalErrorIn
        (
            "Foam::RBCCollection::RBCZoneID(const label RBCIdx) const"
        ) << "Nonexisting cell/point/face zones for RBC " << RBCIdx 
          << "." << nl
          << "Cell zones:  " << RBCMesh_.cellZones() .names() << nl
          << "Point zones: " << RBCMesh_.pointZones().names() << nl
          << "Face zones:  " << RBCMesh_.faceZones() .names() << nl
          << abort(FatalError);
    }

    return -1;
}


bool
RBCCollection::RBCHasNonEmptyZone(const label RBCIdx) const
{
    const word czName = RBCMeshFactory:: cellZoneName(RBCIdx);
    const word pzName = RBCMeshFactory::pointZoneName(RBCIdx);
    const word fzName = RBCMeshFactory:: faceZoneName(RBCIdx);

    const label czID = RBCMesh_.cellZones() .findZoneID(czName);
    const label pzID = RBCMesh_.pointZones().findZoneID(pzName);
    const label fzID = RBCMesh_.faceZones() .findZoneID(fzName);

    if(czID >= 0 && pzID >= 0 && fzID >= 0 && czID == pzID && pzID == fzID)
    {
        return RBCMesh_.hasNonEmptyZone(czID);
    }
    return false;
}


void
RBCCollection::setConstantPO2(const scalar PO2_RBC)
{
    PO2_ = PO2_RBC;
    Hb_ = dissociationCurve_.equilibriumS(PO2_);
}


void
RBCCollection::setConstantPO2(const dimensionedScalar PO2_RBC)
{
    setConstantPO2(PO2_RBC.value());
}


void
RBCCollection::setConstantPO2(const label RBCIdx, const scalar PO2_RBC)
{
    const scalar eqHb = dissociationCurve_.equilibriumS(PO2_RBC);

    RBCMesh_.setFieldOnZone(PO2_, PO2_RBC, RBCZoneID(RBCIdx)); 
    RBCMesh_.setFieldOnZone(Hb_,     eqHb, RBCZoneID(RBCIdx)); 
}


void
RBCCollection::setConstantPO2
(
    const label RBCIdx, 
    const dimensionedScalar PO2_RBC
)
{
    setConstantPO2(RBCIdx, PO2_RBC.value());
}


void
RBCCollection::setEquilibriumPO2
(
    const label RBCIdx 
)
{
    scalarField zoneHb (RBCMesh_.getFieldOnZone(Hb_, RBCZoneID(RBCIdx)));
    scalarField zonePO2(dissociationCurve_.equilibriumPO2(zoneHb));
    RBCMesh_.setFieldOnZone
    (
        PO2_,
        zonePO2,
        RBCZoneID(RBCIdx)
    );

    RBCMesh_.setBoundaryFieldToPatchInternalFieldOnZone(PO2_, RBCZoneID(RBCIdx));
}


void
RBCCollection::setEquilibriumHb()
{
    forAll(RBCs_.toc(), i)
    {
        const label RBCIdx = RBCs_.toc()[i];
        if (RBCHasNonEmptyZone(RBCIdx))
        {
            setEquilibriumHb(RBCIdx);
        }
    }
}


void
RBCCollection::setEquilibriumHb
(
    const label RBCIdx 
)
{
    scalarField zonePO2(RBCMesh_.getFieldOnZone(PO2_, RBCZoneID(RBCIdx)));
    scalarField zoneHb(dissociationCurve_.equilibriumS(zonePO2));
    RBCMesh_.setFieldOnZone
    (
        Hb_,
        zoneHb,
        RBCZoneID(RBCIdx)
    );

    RBCMesh_.setBoundaryFieldToPatchInternalFieldOnZone(Hb_, RBCZoneID(RBCIdx));
}


void
RBCCollection::insert
(
    const RBC& cRBC
)
{
    if (debug)
    {
        Info<< "Foam::RBCCollection::insert(const label, "
               "const deformableBodyGeometricState&, "
               "const word&, const scalar):"
            << " Inserting RBC " << cRBC.index() << "." << endl;
    }

    RBCs_.insert(cRBC.index(), autoPtr<RBC>(new RBC(cRBC)));

    addMeshIfInDomain(cRBC);
}


void
RBCCollection::remove(const label RBCIdx)
{
    if (RBCs_.find(RBCIdx) == RBCs_.end())
    {
        FatalErrorIn
        (
            "Foam::RBCCollection::remove(const label)"
        ) << "RBC index " << RBCIdx << " not found."
          << abort(FatalError);
    }

    if (debug)
    {
        Info<< "Foam::RBCCollection::remove(const label): "
            << "Removing RBC " << RBCIdx << endl;
    }

    // if data are present, remove RBC zones
    if (RBCHasNonEmptyZone(RBCIdx))
    {
        removeMesh(RBCs_[RBCIdx]);
    }

    // delete RBC object
    RBCs_.erase(RBCIdx);
}


void
RBCCollection::updateMeshes()
{
    // apply changes to the mesh
    RBCMesh_.update();

    // move the added mesh zones to their positions
    moveMeshes();

    // set initial values on the new mesh zones
    initializeAddedFields();

    // clear the list of added mesh zones
    addedList_.clear();
}


void
RBCCollection::prepareMotion
(
    const label RBCIdx,
    const deformableBodyGeometricState& geometricState
)
{
    RBC& cRBC = (*this)[RBCIdx];
    // store the old RBC geometric state
    autoPtr<deformableBodyGeometricState> oldGeometricStatePtr = cRBC.geometricState().clone();
    const deformableBodyGeometricState& oldGeometricState = oldGeometricStatePtr();
    // replace the geometric state
    cRBC.setGeometricState(geometricState);

    if (Pstream::parRun())
    {
        const label currentProcI = deformableBodyProc(oldGeometricState);
        const label nextProcI    = deformableBodyProc(   geometricState);

        bool nowInDomain  = (currentProcI >= 0);
        bool nextInDomain = (   nextProcI >= 0);
        // List that contain all the indices of RBCs that are to be sent
        // to/received from the current processor
        labelListList sendList(Pstream::nProcs());
        labelListList recvList(Pstream::nProcs());

        if (debug & 2)
        {
            Info<< "Foam::RBCCollection::prepareMotion(const label, "
                << "const point&, const vector&): "
                << "RBC " << RBCIdx << ": currentProcI = " << currentProcI
                << ", nextProcI = " << nextProcI << endl;
        }

        // If the RBC goes from current processor to outside the domain
        if ((currentProcI == Pstream::myProcNo()) && (!nextInDomain))
        {
            removeMesh(cRBC);
        }
        // If the RBC goes from outside the domain to current processor
        else if (!nowInDomain && nextProcI == Pstream::myProcNo())
        {
            addMesh(cRBC);
            addedList_.append(RBCIdx);
        }
        // If the RBC goes from current processor to other processor
        else if 
        ( 
            (currentProcI == Pstream::myProcNo()) &&
            (nextProcI    != Pstream::myProcNo()) &&
            (nextInDomain)
        )
        {
            Pout<< "Foam::RBCCollection::prepareMotion(const label, "
                << "const point&, const vector&): "
                << "Going to send RBC " << RBCIdx << endl;
            removeMesh(cRBC);
            sendList_[nextProcI].append(RBCIdx);
        }
        // If the RBC goes from other processor to current processor
        else if 
        ( 
            (currentProcI != Pstream::myProcNo()) &&
            (nowInDomain)                         &&
            (nextProcI    == Pstream::myProcNo())
        )
        {
            Pout<< "Foam::RBCCollection::prepareMotion(const label, "
                << "const point&, const vector&): "
                << "Going to receive RBC " << RBCIdx << endl;
            addMesh(cRBC);
            recvList_[currentProcI].append(RBCIdx);
        }
    }
    else
    {
        bool nowInDomain  = deformableBodyInDomain(oldGeometricState);
        bool nextInDomain = deformableBodyInDomain(   geometricState);
        // if RBC leaves domain
        if (nowInDomain && (!nextInDomain))
        {
            removeMesh(cRBC);
        }
        // if RBC enters domain
        else if ((!nowInDomain) && nextInDomain)
        {
            addMesh(cRBC);
            addedList_.append(RBCIdx);
        }
    }
}


void
RBCCollection::applyMotion()
{
    sendFields();

    // apply changes to the mesh, map fields, move meshes and initialize fields
    updateMeshes();
    
    receiveFields();
    setOldReceivedFields();

    // set PO2_RBC to equilibrium and set boundary values of Hb and PO2.
    forAll(recvList_, procI)
    {
        labelList RBCList = recvList_[procI];
        forAll(RBCList, i)
        {
            const label RBCIdx = RBCList[i];
            RBCMesh_.setBoundaryFieldToPatchInternalFieldOnZone(Hb_, RBCZoneID(RBCIdx));
            setEquilibriumPO2(RBCIdx);
            setRBCVolumeFraction(RBCIdx);
        }
    }
    
    // clear lists for communicated RBCs
    forAll(sendList_, i)
    {
        sendList_[i].clear();
        recvList_[i].clear();
    }
}


scalar 
RBCCollection::PO2FromFieldNameAndValue(const word& fieldName, const scalar fieldValue)
{
    if (fieldName == "PO2")
    {
        return fieldValue;
    }
    else if (fieldName == "Hb")
    {
        return dissociationCurve_.equilibriumPO2(fieldValue);
    }
    else
    {
        FatalErrorIn
        (
            "Foam::RBCCollection::PO2FromFieldNameAndValue(const word&, const scalar):"
        ) << "Invalid field name " << fieldName << ". " << nl
          << "Supported field names are PO2 and Hb."
          << abort(FatalError);
    }
    return 0.0;
}


bool
RBCCollection::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


bool
RBCCollection::write() const
{
    return regIOobject::write();
}


bool
RBCCollection::read()
{
    return regIOobject::read();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

const RBC&
RBCCollection::operator[](const label i) const
{
    if (!RBCs_.found(i))
    {
        FatalErrorIn
        (
            "Foam::RBCCollection::operator[](const label) const"
        ) << "There is no RBC with index " << i << "." << nl
          << abort(FatalError);
    }
    return RBCs_[i]();
}


RBC&
RBCCollection::operator[](const label i)
{
    if (!RBCs_.found(i))
    {
        FatalErrorIn
        (
            "Foam::RBCCollection::operator[](const label) const"
        ) << "There is no RBC with index " << i << "." << nl
          << abort(FatalError);
    }
    return RBCs_[i]();
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const RBCCollection& RBCColl
)
{
    forAllConstIter(Map<autoPtr<RBC> >, RBCColl.RBCs_, iter)
    {
        const RBC& cRBC = *iter;
        os  << "RBC" << cRBC.index() << nl
            << token::BEGIN_BLOCK << incrIndent << nl 
            << cRBC
            << decrIndent << indent << token::END_BLOCK << endl;
    }
    return os;
}


// ************************************************************************* //
