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

\*---------------------------------------------------------------------------*/

#include "implicitSourceControl.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(implicitSourceControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::implicitSourceControl::read()
{
    solutionControl::read(true);

    // Read solution controls
    const dictionary& implicitSourceDict = dict();
    nIterMax_         = implicitSourceDict.lookupOrDefault<label>("nIterationsMax", 1000);

    // sanity check
    if (nIterMax_ < 1)
    {
        FatalErrorIn
        (
            "implicitSourceControl::read()"
        )   << " nIterationsMax should be greater than zero." << nl
            << abort(FatalError);
    }
}


bool Foam::implicitSourceControl::criteriaSatisfied()
{
    // no checks on first iteration - nothing has been calculated yet
    if ((corr_ == 1) || residualControl_.empty() || finalIter())
    {
        return false;
    }

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    const dictionary& solverDict = mesh_.solverPerformanceDict();
    forAllConstIter(dictionary, solverDict, iter)
    {
        const word& variableName = iter().keyword();
        const label fieldI = applyToField(variableName);
        if (fieldI != -1)
        {
            const List<solverPerformance> sp(iter().stream());
            const scalar residual = sp.last().initialResidual();

            checked = true;

            const bool absCheck = residual < residualControl_[fieldI].absTol;
            bool relCheck = false;

            achieved = achieved && (absCheck || relCheck);

            if (debug)
            {
                Info<< algorithmName_ << " loop:" << endl;

                Info<< "    " << variableName
                    << " Implicit source term iteration " << corr_
                    << ": ini res = "
                    << residualControl_[fieldI].initialResidual
                    << ", abs tol = " << residual
                    << " (" << residualControl_[fieldI].absTol << ")"
                    << endl;
            }
        }
    }

    return checked && achieved;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::implicitSourceControl::implicitSourceControl(fvMesh& mesh)
:
    solutionControl(mesh, "implicitSource"),
    nIterMax_(1)
{
    read();

    if (nIterMax_ > 1)
    {
        Info<< nl;
        if (residualControl_.empty())
        {
            Info<< algorithmName_ << ": no residual control data found. "
                << "Calculations will employ " << nIterMax_
                << " corrector loops" << nl << endl;
        }
        else
        {
            Info<< algorithmName_ << ": max iterations = " << nIterMax_
                << endl;
            forAll(residualControl_, i)
            {
                Info<< "    field " << residualControl_[i].name << token::TAB
                    // << ": relTol " << residualControl_[i].relTol
                    << ": tolerance " << residualControl_[i].absTol
                    << nl;
            }
            Info<< endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::implicitSourceControl::~implicitSourceControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::implicitSourceControl::loop()
{
    read();

    corr_++;

    if (debug)
    {
        Info<< algorithmName_ << " loop: corr = " << corr_ << endl;
    }

    if (corr_ == nIterMax_ + 1)
    {
        if ((!residualControl_.empty()) && (nIterMax_ != 1))
        {
            Info<< algorithmName_ << ": not converged within "
                << nIterMax_ << " iterations" << endl;
        }

        corr_ = 0;
        mesh_.data::remove("finalIteration");
        return false;
    }

    bool completed = false;
    if (criteriaSatisfied())
    {
        Info<< algorithmName_ << ": converged in " << corr_ - 1
            << " iterations" << endl;

        mesh_.data::remove("finalIteration");
        corr_ = 0;

        completed = true;
    }
    else
    {
        if (finalIter())
        {
            mesh_.data::add("finalIteration", true);
        }

        if (corr_ <= nIterMax_)
        {
            if (nIterMax_ != 1)
            {
                Info<< algorithmName_ << ": iteration " << corr_ << endl;
                storePrevIterFields();
            }

            completed = false;
        }
    }

    return !completed;
}

bool Foam::implicitSourceControl::finalIter() const
{
    return (corr_ == nIterMax_);
}

// ************************************************************************* //

