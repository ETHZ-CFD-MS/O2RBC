/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::implicitSourceControl

Description
    Control class for the linearization of implicit source terms. The iteration
    stops when the difference between the current field and thr previous
    iteration is less than a threshold, or when the maximum number of
    iterations is reached.

\*---------------------------------------------------------------------------*/

#ifndef implicitSourceControl_H
#define implicitSourceControl_H

#include "solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class implicitSourceControl Declaration
\*---------------------------------------------------------------------------*/

class implicitSourceControl
:
    public solutionControl
{
protected:

    // Protected data

        // Solution controls

            //- Maximum number of iterations
            label nIterMax_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual void read();

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();

        //- Disallow default bitwise copy construct
        implicitSourceControl(const implicitSourceControl&);

        //- Disallow default bitwise assignment
        void operator=(const implicitSourceControl&);


public:


    // Static Data Members

        //- Run-time type information
        TypeName("implicitSourceControl");


    // Constructors

        //- Construct from mesh
        implicitSourceControl(fvMesh& mesh);


    //- Destructor
    virtual ~implicitSourceControl();


    // Member Functions

        // Access

            //- Maximum number of iterations
            label nIterMax() const
            {
                return nIterMax_;
            }


        // Solution control

            //- Loop control
            virtual bool loop();

            //- Helper function to identify final iteration
            bool finalIter() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

