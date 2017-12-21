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

Class
    Foam::RBCCollection

Description
    Container for RBCs.

SourceFiles
    RBCCollection.C

\*---------------------------------------------------------------------------*/

#ifndef RBCCollection_H
#define RBCCollection_H

#include "RBC.H"
#include "disconnectedZoneMesh.H"
#include "cartesianProcMeshInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class OFstream;
class dissociationCurve;


// Forward declaration of friend functions and operators
class RBCCollection;
Ostream& operator<<(Ostream&, const RBCCollection&);


/*---------------------------------------------------------------------------*\
                         Class RBCCollection Declaration
\*---------------------------------------------------------------------------*/

class RBCCollection
: 
    public regIOobject
{
    // Private data

        //- RBC mesh
        disconnectedZoneMesh RBCMesh_;
    
        //- Hb field
        volScalarField Hb_;

        //- PO2 field
        volScalarField PO2_;
        
        //- Indicator field of RBC in RBC mesh
        volScalarField in_RBC_;

        //- Reference to dissociation curve object.
        const dissociationCurve& dissociationCurve_;

        //- Map of RBC objects
        Map<autoPtr<RBC> > RBCs_;

        //- Information on processor meshes
        const cartesianProcMeshInfo procMeshInfo_;

        //- Lists of RBCs to send to processors
        labelListList sendList_;

        //- Lists of RBCs to receive from processors
        labelListList recvList_;

        //- List of RBCs that are added "from scratch"
        labelList addedList_;

        //- Hash table of blue print meshes
        mutable HashTable<autoPtr<polyMesh>, word, Hash<word> > bluePrintMeshTable_;

        //- Hash table ofbase geometric states
        mutable HashTable
        <
            autoPtr<deformableBodyGeometricState>, 
            word, 
            Hash<word> 
        > geomStateTable_;


    // Private Member Functions

        //- Add a RBC to the RBC mesh
        void addMesh(const RBC&);

        //- Add a RBC to the RBC mesh if it is inside the domain
        void addMeshIfInDomain(const RBC&);

        //- Add a RBC to the RBC mesh
        void removeMesh(const RBC&);

        //- Move the mesh zones to their associated RBC position
        void moveMeshes();

        //- Initialize fields on the added zones
        void initializeAddedFields();

        //- Set the volume fraction field on a zone
        void setRBCVolumeFraction(const label RBCIdx);

        //- Set old time values on a given RBC zone for fields for which an 
        //  equation is solved
        void setFieldOldTime(const label RBCIdx);

        //- Set old time values on a given RBC zone for the mesh volume
        void setVolumeOldTime(const label RBCIdx);

        //- Send fields for RBC zones that leave the current processor
        void sendFields();

        //- Receive fields for RBC zones that enter the current processor
        void receiveFields();

        //- Set the old values of the received fields.
        void setOldReceivedFields();

        //- Return whether the RBC overlaps the computational domain
        bool RBCInDomain(const label RBCIdx) const;

        //- Return whether the given body overlaps the computational domain
        bool deformableBodyInDomain(const deformableBodyGeometricState&) const;

        //- Return the processor number that a RBC belongs to geometrically
        label geometricRBCProc(const label RBCIdx) const;

        //- Return the processor number that a deformable body belongs to geometrically
        label deformableBodyProc(const deformableBodyGeometricState&) const;

        //- Add a blue print mesh to the hash table with sample mesh names as
        //  key
        void addMeshToTable(const word& sampleMeshName, autoPtr<polyMesh>);

        //- Return the mesh associated to a sample mesh name
        const polyMesh& bluePrintMesh(const word& sampleMeshName) const;

        //- Return the geometric state associated to a sample mesh name
        // autoPtr<deformableBodyGeometricState> geometricState
        const deformableBodyGeometricState& geometricState
        (
            const word& sampleMeshName
        ) const;

        //- Disallow default bitwise copy construct
        RBCCollection(const RBCCollection&);

        //- Disallow default bitwise assignment
        void operator=(const RBCCollection&);


public:

    //- Run-time type information
    TypeName("dictionary"); // writing and reading do not work if RBCCollection


    // Constructors

        //- Construct using mesh object and dissociationCurve.
        RBCCollection
        (
            const fvMesh& mesh, 
            const dissociationCurve&
        );

    //- Destructor
    ~RBCCollection();


    // Member Functions

        // Access
            
            //- Return number of RBCs
            label size() const
            {
                return RBCs_.size();
            }

            labelList RBCIndices() const
            {
                return RBCs_.toc();
            }

            //- Constant access to the RBC mesh
            const disconnectedZoneMesh&
            RBCMesh() const
            {
                return RBCMesh_;
            }

            //- Constant access to PO2
            const volScalarField& 
            PO2() const
            {
                return PO2_;
            }

            //- Non-constant access to PO2
            volScalarField& 
            PO2()
            {
                return PO2_;
            }
            
            //- Constant access to Hb
            const volScalarField& 
            Hb() const
            {
                return Hb_;
            }

            //- Non-constant access to Hb
            volScalarField& 
            Hb()
            {
                return Hb_;
            }
            
            //- Constant access to in_RBC
            const volScalarField& 
            in_RBC() const
            {
                return in_RBC_;
            }

            //- Non-constant access to a field with given name
            //  Currently, it is allowed for Hb and PO2 only.
            volScalarField&
            getField(const word& fieldName); 
           
            //- Constant access to the dissociation curve
            const dissociationCurve&
            getDissociationCurve() const
            {
                return dissociationCurve_;
            }

            //- Return the zoneID for a given RBC.
            //  Return -1 if the given RBC has no corresponding zone.
            label RBCZoneID(const label RBCIdx) const;

            //- Return whether a RBC has associated nonempy mesh zones
            bool RBCHasNonEmptyZone(const label RBCIdx) const;

        // Setters

            //- Set a constant PO2 value in all RBCs.
            //  Set hemoglobin in equilibrium with PO2_RBC.
            void setConstantPO2(const scalar PO2_RBC);

            //- Set a constant PO2 value in the all RBCs.
            //  Set hemoglobin in equilibrium with PO2_RBC.
            void setConstantPO2(const dimensionedScalar PO2_RBC);

            //- Set a constant PO2 value in the i-th RBC.
            //  Set hemoglobin in equilibrium with PO2_RBC.
            void setConstantPO2(const label, const scalar PO2_RBC);

            //- Set a constant PO2 value in the i-th RBC and in the plasma.
            //  Set hemoglobin in equilibrium with PO2_RBC.
            void setConstantPO2(const label, const dimensionedScalar PO2_RBC);

            //- Set PO2 in equilibrium with the hemoblogin saturation in a given
            //  RBC.
            void setEquilibriumPO2(const label);

            //- Set hemoglobin saturation in equilibrium with PO2 in all RBCs
            void setEquilibriumHb();

            //- Set hemoglobin saturation in equilibrium with PO2 in a given RBC
            void setEquilibriumHb(const label);


        // Edit

            //- Insert a RBC into the collection with an inlet field name 
            //  and an inlet value.
            //  Prepares the region creation, but does not 
            //  perform mesh creation and field data initialization.
            void insert(const RBC&);

            //- Remove a RBC
            void remove(const label);

            //- Creates meshes. Needs to be called after RBC insertion or
            //  removal.
            void updateMeshes();

            //- Prepare the motion of a RBC. This does not change the mesh
            //  topology, nor does parallel communication
            void prepareMotion
            (
                const label RBCIdx, 
                const deformableBodyGeometricState& geometricState
            );

            //- Apply motion by changing the mesh and doing parallel
            //  communication.
            void applyMotion();

        // Miscellaneous

            //- Return PO2 from a field name and a value
            scalar PO2FromFieldNameAndValue(const word&, const scalar);

        // Write

            //- Write function used by regIOobject for writing
            bool writeData(Ostream& os) const;

            //- Write
            bool write() const;

        // Read

            bool read();

        // Operators
            const RBC& operator[](const label) const;

            RBC& operator[](const label);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const RBCCollection&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
