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
    Foam::tubeVisitor

Description
    Visitor abstract class for the class circularTube.

SourceFiles
    tubeVisitor.C

\*---------------------------------------------------------------------------*/

#ifndef tubeVisitor_H
#define tubeVisitor_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polygonalTube;

/*---------------------------------------------------------------------------*\
                         Class tubeVisitor Declaration
\*---------------------------------------------------------------------------*/

class tubeVisitor
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        tubeVisitor(const tubeVisitor&);

        //- Disallow default bitwise assignment
        void operator=(const tubeVisitor&);


public:

    // Constructors

        //- Construct null
        tubeVisitor();

    // Destructor

        virtual ~tubeVisitor() {}

    // Member Functions

        virtual void visitPolygonalTube(const polygonalTube&) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
