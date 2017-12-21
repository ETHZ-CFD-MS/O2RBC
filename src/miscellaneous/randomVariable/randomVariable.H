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
    Foam::randomVariable

Description
    Contains parameters and sampling for various families of random variables. 
    Currently supports constant, uniform and normal random variables. A list
    of deterministic values taken in a cyclic manner is also provided.

    The class Random from the OpenFOAM library is not used since the random
    number generation in parallel was sometimes not working with that class.

SourceFiles
    randomVariable.H
    randomVariable.C

\*---------------------------------------------------------------------------*/

#ifndef randomVariable_H
#define randomVariable_H

#include "fvCFD.H"

#include <cstdlib>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class randomVariable Declaration
\*---------------------------------------------------------------------------*/

class randomVariable
{

public:
    enum randomType
    {
        CONSTANT,
        CYCLIC_LIST,
        DISCRETE,
        UNIFORM, 
        GAUSSIAN,
        LOG_NORMAL_SPACING
    };

private:
    // Private data

        //- Type of random variable
        randomType type_;

        //- Weights (for discrete random variable)
        scalarList weights_;

        //- Values (for cyclic list of discrete random variable)
        scalarList values_;

        //- Lower bound for the values (for sanity checks)
        scalar lowerBound_;

        //- Upper bound for the values (for sanity checks)
        scalar upperBound_;

        //- Mean value (for Gaussian).
        scalar meanValue_;

        //- Standard deviation (for Gaussian)
        scalar standardDev_;

        //- Maximum value (for uniform and truncated Gaussian)
        scalar maxValue_;

        //- Minimum value (for uniform and truncated Gaussian)
        scalar minValue_;

        //- mu parameter for log-normal spacing
        scalar mu_;

        //- sigma parameter for log-normal spacing
        scalar sigma_;

        //- counter for the number of calls to generate()
        mutable label generateCounter_;


    // Private Member Functions
    
        //- Read attributes from dictionary
        bool read(const dictionary& dict);

        //- Initialize random number generator
        void initializeRandomGenerator() const;

        //- Generate a seed
        int seed() const;

        //- Uniform random number between 0 and 1
        scalar scalar01() const;

        //- Normal random number with mean 0 and variance 1
        scalar GaussNormal() const;

        //- Perform sanity checks on minValue_ and maxValue_
        void checkMinMax() const;

        //- Perform sanity checks on mean value
        void checkMean() const;

        //- Perform sanity checks on standard deviation
        void checkStdDev() const;

        //- Disallow default bitwise assignment
        void operator=(const randomVariable&);


public:

    // Constructors

        //- Default constructor
        randomVariable();

        //- Construct from name and dictionary
        randomVariable
        (
            const dictionary& dict
        );

        //- Construct from name, dictionary and bounds
        randomVariable
        (
            const dictionary& dict,
            const scalar lowerBound,
            const scalar upperBound
        );

        //- Construct from Istream
        randomVariable
        (
            Istream& is
        );

        //- Copy constructor
        randomVariable(const randomVariable&);


    //- Destructor
    ~randomVariable();


    // Member Functions

        // Access

            randomType type() const
            {
                return type_;
            }

            scalar meanValue() const;

        // Evaluation

            //- Generate a random sample
            scalar generate() const;

        // Write
        
            void write(Ostream& os) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, randomVariable&);
        friend Ostream& operator<<(Ostream&, const randomVariable&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
