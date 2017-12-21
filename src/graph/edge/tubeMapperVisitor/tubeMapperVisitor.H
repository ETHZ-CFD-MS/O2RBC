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

Class
    Foam::tubeMapperVisitor

Description
    Class for the construction of a polyMesh object around a tube of any type.

    This code assumes that the cylinder axis is along the x-direction.

SourceFiles
    tubeMapperVisitor.C

\*---------------------------------------------------------------------------*/

#ifndef tubeMapperVisitor_H
#define tubeMapperVisitor_H

#include "tubeVisitor.H"

#include "fvCFD.H"

#include "label.H"
#include "scalar.H"
#include "List.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "boundBox.H"

#include "polygonalTube.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;

// Forward declaration of friend functions and operators
class tubeMapperVisitor;
Ostream& operator<<(Ostream&, const tubeMapperVisitor&);


/*---------------------------------------------------------------------------*\
                         Class tubeMapperVisitor Declaration
\*---------------------------------------------------------------------------*/

class tubeMapperVisitor
:
    public tubeVisitor
{
    // Private data

        //- Cylinder mesh used as a base mesh.
        fvMesh cylinderMesh_;

        //- Transformed points for the mapped mesh creation.
        pointField transformedPoints_;

        //- Cylinder axis points
        Pair<point> cylinderAxisPoints_;

        //- Cylinder radius
        scalar cylinderRadius_;

        //- List of rotation tensors that map segments to the cylinder axis
        List<tensor> rotationToAxisList_;

    // Private Member Functions

        //- Compute the geometric data of the cylinder mesh
        void computeCylinderData();

        //- Compute all the geometric transformation needed the mapping
        void computeGeometricTransformations(const polygonalTube&);

        //- Disallow default bitwise copy construct
        tubeMapperVisitor(const tubeMapperVisitor&);

        //- Disallow default bitwise assignment
        void operator=(const tubeMapperVisitor&);


public:

    // Runtime type information
    ClassName("tubeMapperVisitor");

    // Static data members

        //- Static data staticData
        static const word tubeMeshPrefix_;

    // Static member functions

        //- Returns the tube mesh name for a given index
        static word tubeMeshName(const label);

    // Constructors

        //- Construct from IOobject
        tubeMapperVisitor(const IOobject&);


    //- Destructor
    ~tubeMapperVisitor();


    // Member Functions

        // Access
            
            //- Constant access to the cylinder mesh
            const fvMesh& cylinderMesh() const
            {
                return cylinderMesh_;
            }

            //- Access to the mapped mesh.
            //  The argument is used to create the mesh name.
            autoPtr<polyMesh> mappedMesh(const label) const;

            //- Return a field with the normalized axis-parallel of the cylinder
            //  cell centers.
            autoPtr<volScalarField> normalizedCylinderAxisCoord() const;

        // Visitor method
        
            //- Create a mapped mesh given a polygonal tube
            virtual void visitPolygonalTube(const polygonalTube&);

        // Geometric constructions
            
            //- Transform a point to its position in the tube mesh
            point transformPoint
            (
                const point&,
                const polygonalTube&,
                const label segmentI
            ) const;

            //- Construct the rotation tensor that maps the segment axis to the
            //  cylinder axis.
            tensor rotationTensorFromSegmentToAxis
            (
                const polygonalTube&,
                const label segmentI
            ) const;

            //- Construct the shear tensor for given ellipse axes and the
            //  corresponding rotation that maps it to the cylinder axis. 
            tensor shearTensorOnCylinderAxis
            (
                const ellipseAxes&, 
                const tensor&
            ) const;

            //- Returns the indices of the cylinder point that will be mapepd
            //  to the current segment.
            labelList cylinderPointIDsOnSegment
            (
                const polygonalTube& tube,
                const label segmentI
            ) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const tubeMapperVisitor&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
