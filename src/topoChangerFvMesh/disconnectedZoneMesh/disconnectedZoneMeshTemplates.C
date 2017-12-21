/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::Field<Type>
Foam::disconnectedZoneMesh::getFieldOnZone
(
    const Field<Type>& f,
    const label zoneI
) const
{
    const Field<Type> zoneField(f, cellZones()[zoneI]);
    return zoneField;
}


template<class Type>
Foam::Field<Type>
Foam::disconnectedZoneMesh::getFieldOnZone
(
    const GeometricField<Type,fvPatchField,volMesh>& f, 
    const label zoneI
) const
{
    return getFieldOnZone(f.field(), zoneI);
}


template<class Type>
void Foam::disconnectedZoneMesh::setFieldOnZone
(
    Field<Type>& f, 
    const Field<Type>& zoneField,
    const label zoneI
) const
{
    const cellZone& cz = cellZones()[zoneI];
    if (zoneField.size() != cz.size())
    {
        FatalErrorIn
        (
            "Foam::disconnectedZoneMesh::setFieldOnZone(Field<Type>&, "
            "const Field<Type>&, const label) const"
        ) << "The field zoneField and the zone " << zoneI
          << " do not have the same size."
          << abort(FatalError);
    }
    forAll(cz, i)
    {
        f[cz[i]] = zoneField[i];
    }
}


template<class Type>
void Foam::disconnectedZoneMesh::setFieldOnZone
(
    Field<Type>& f, 
    const Type value,
    const label zoneI
) const
{
    Field<Type> zoneField(cellZones()[zoneI].size(), value);
    setFieldOnZone(f, zoneField, zoneI);
}


template<class Type>
void Foam::disconnectedZoneMesh::setFieldOnZone
(
    GeometricField<Type,fvPatchField,volMesh>& f, 
    const Type value,
    const label zoneI
) const
{
    setFieldOnZone(f.primitiveFieldRef(), value, zoneI);
    setFieldOnZoneBoundary(f, value, zoneI);
}


template<class Type>
void Foam::disconnectedZoneMesh::setFieldOnZoneBoundary
(
    GeometricField<Type,fvPatchField,volMesh>& f, 
    const Type value,
    const label zoneI
) const
{
    const labelPairList zoneBoundaryIndices = zoneBoundaryFaces(zoneI);
    forAll(zoneBoundaryIndices, i)
    {
        label patchI = zoneBoundaryIndices[i].first();
        label faceI  = zoneBoundaryIndices[i].second();
        f.boundaryFieldRef()[patchI][faceI] = value;
    }
}


template<class Type>
void Foam::disconnectedZoneMesh::setBoundaryFieldToPatchInternalFieldOnZone
(
    GeometricField<Type,fvPatchField,volMesh>& f, 
    const label zoneI
) const
{
    const labelPairList zoneBoundaryIndices = zoneBoundaryFaces(zoneI);
    forAll(zoneBoundaryIndices, i)
    {
        label patchI = zoneBoundaryIndices[i].first();
        label faceI  = zoneBoundaryIndices[i].second();
        label cellI  = boundary()[patchI].faceCells()[faceI];
        f.boundaryFieldRef()[patchI][faceI] = f[cellI];
    }
}

