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

#include "RBCProperty.H"
#include "graphCoordinateDirected.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * *  Private member functions * * * * * * * * * * * * //

template<class T>
void Foam::RBCProperty<T>::checkType() const
{
    if (type_ != "constant" &&
        type_ != "edgeFunction")
    {
        FatalErrorIn("Foam::RBCProperty<T>::checkType()")
            << "Unknown type of RBCProperty " << type_ << "." << nl
            << "Valid types are ""constant"" and ""edgeFunction""." << nl
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::RBCProperty<T>::RBCProperty
(
    const dictionary& dict
)
:   type_(dict.lookup("type"))
{
    checkType();
    if (type_ == "constant")
    {
        constantValue_ = T(dict.lookup("value"));
    }
    else if (type_ == "edgeFunction")
    {
        edgeValues_ = dict.lookup("values");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
List<T> RBCProperty<T>::values() const
{
    if (type_ == "constant")
    {
        return List<T>(1, constantValue_);
    }
    else if (type_ == "edgeFunction")
    {
        List<T> valueList();
        forAllConstIter(typename Map<T>, edgeValues_, iter)
        {
            valueList.append(*iter);
        }
        return valueList;
    }
}


template<class T>
T RBCProperty<T>::value(const label RBCIdx, const label edgeI) const
{
    if (type_ == "constant")
    {
        return constantValue_;
    }
    else if (type_ == "edgeFunction")
    {
        return edgeValues_[edgeI];
    }
    return constantValue_;  // to avoid compiler warning
}


template<class T>
void RBCProperty<T>::write(Ostream& os) const
{
    os.writeKeyword("type")   << type_ << token::END_STATEMENT << nl;
    if (type_ == "constant")
    {
        os.writeKeyword("value") << constantValue_ << token::END_STATEMENT << nl;
    }
    else if (type_ == "edgeFunction")
    {
        os.writeKeyword("values") << edgeValues_ << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

template<class T>
Foam::Ostream&
Foam::operator<<(Ostream& os, const RBCProperty<T>& obj)
{
    obj.write(os);

    return os;
}


// ************************************************************************* //
