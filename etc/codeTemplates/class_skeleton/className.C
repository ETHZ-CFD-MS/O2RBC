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

#include "className.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

className::className()
:
    attr1_(),
    attr2_()
{}

className::className(scalar attr1, scalar attr2)
:
    attr1_(attr1),
    attr2_(attr2)
{}

className::className
(
    const className& obj
)
:
    attr1_(obj.attr1_),
    attr2_(obj.attr2_)
{}

void className::write(Ostream& os) const
{
    os << "attr1 = " << attr1_ << endl;
    os << "attr2 = " << attr2_ << endl;
}




// ************************************************************************* //
