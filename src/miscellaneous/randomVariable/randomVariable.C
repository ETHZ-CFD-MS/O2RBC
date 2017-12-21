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

#include "randomVariable.H"

// For parallel debugging
#include "unistd.h"
#include "stdio.h"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool 
randomVariable::read(const dictionary& dict)
{
    word typeWord = dict.lookup("type");
    if (typeWord == "constant")
    {
        type_         = CONSTANT;
        meanValue_    = readScalar(dict.lookup("value"));
    }
    else if (typeWord == "cyclicList")
    {
        type_         = CYCLIC_LIST;
        values_       = dict.lookup("values");
    }
    else if (typeWord == "discrete")
    {
        type_         = DISCRETE;
        weights_      = dict.lookup("weights");
        values_       = dict.lookup("values");
        if (weights_.size() != values_.size())
        {
            FatalErrorIn
            (
                "randomVariable::read(const dictionary& dict)"
            )   << "  the lists of weights and values do not have the same length."
                << abort(FatalError);
        }
        scalar weightSum = Foam::sum(weights_);
        forAll(weights_, i)
        {
            weights_[i] /= weightSum;
        }
    }
    else if (typeWord == "uniform")
    {
        type_         = UNIFORM;
        minValue_     = readScalar(dict.lookup("minValue"));
        maxValue_     = readScalar(dict.lookup("maxValue"));
        meanValue_    = 0.5*(minValue_ + maxValue_);
        checkMinMax();
    }
    else if (typeWord == "gaussian")
    {
        type_         = GAUSSIAN;
        meanValue_    = readScalar(dict.lookup("mean"));
        standardDev_  = readScalar(dict.lookup("stdDev"));
        minValue_     = dict.lookupOrDefault<scalar>("minValue", -GREAT);
        maxValue_     = dict.lookupOrDefault<scalar>("maxValue",  GREAT);
        checkMinMax();
        checkMean();
        checkStdDev();

        // do sanity check on mean and standard deviation
        if (meanValue_   <= lowerBound_ || 
            meanValue_   >= upperBound_ || 
            standardDev_ <= 0.0)
        {
            FatalErrorIn
            (
                "randomVariable::read(const dictionary& dict)"
            )   << "  invalid value for Mean or StdDev." << nl
                << "  The values should satisfy lowerBound < mean < upperBound and " << nl
                << "  0 < stdDev" << nl
                << abort(FatalError);
        }
    }
    else if (typeWord == "logNormalSpacing")
    {
        type_         = LOG_NORMAL_SPACING;
        mu_           = readScalar(dict.lookup("mu"));   // mean value of the underlying normal RV
        sigma_        = readScalar(dict.lookup("sigma")); // std. dev. of the underlying normal RV

        // compute approximations to mean value and standard deviation
        // For the mean, compute the median.
        // For the standard deviation, compute a hack.
        // An exact computation would require the evaluation of an integral.
        meanValue_   = 1.0/(1.0 + exp(mu_));
        standardDev_ = 0.5*(1.0/(1.0 + exp(mu_ - sigma_)) + 1.0/(1.0 + exp(mu_ + sigma_)));
        checkStdDev();
    }
    else
    {
        FatalErrorIn
        (
            "randomVariable::read(const dictionary& dict)"
        )   << "  invalid value \"" << typeWord <<  "\" for type." << nl
            << "  valid values are \"constant\", \"cyclicList\", \"discrete\", " << nl
            << "\"uniform\", \"gaussian\" and \"logNormalSpacing\"." << nl
            << abort(FatalError);
    }

    return true;
}


void randomVariable::initializeRandomGenerator() const
{
    srand(seed());
}


int randomVariable::seed() const
{
    return 1 + Pstream::myProcNo();
}


scalar randomVariable::scalar01() const
{
    return static_cast<scalar>(rand())/RAND_MAX;
}


scalar randomVariable::GaussNormal() const
{
    static int iset = 0;
    static scalar gset;
    scalar fac, rsq, v1, v2;

    if (iset == 0)
    {
        do
        {
            v1 = 2.0*scalar01() - 1.0;
            v2 = 2.0*scalar01() - 1.0;
            rsq = v1*v1 + v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);

        fac = sqrt(-2.0 * log(rsq)/rsq);
        gset = v1*fac;
        iset = 1;

        return v2*fac;
    }
    else
    {
        iset = 0;

        return gset;
    }
}


void
randomVariable::checkMinMax() const
{
    if (minValue_ < lowerBound_ || maxValue_ > upperBound_ || minValue_ > maxValue_)
    {
        FatalErrorIn
        (
            "randomVariable::read(const dictionary& dict)"
        )   << "  invalid value for minValue_ or maxValue_." << nl
            << "  The values should satisfy " << nl
            << "  " << lowerBound_ << " < " << minValue_ << " <= " 
            << maxValue_ << " <= " << upperBound_ << nl
            << abort(FatalError);
    }
}

void
randomVariable::checkMean() const
{
    if (meanValue_ <= lowerBound_ || meanValue_ >= upperBound_)
    {
        FatalErrorIn
        (
            "randomVariable::read(const dictionary& dict)"
        )   << "  invalid value for mean." << nl
            << "  The value should satisfy lowerBound < mean < upperBound." << nl
            << abort(FatalError);
    }
}

void
randomVariable::checkStdDev() const
{
    if (standardDev_ <= 0.0)
    {
        FatalErrorIn
        (
            "randomVariable::read(const dictionary& dict)"
        )   << "  invalid value for stdDev." << nl
            << "  The value should satisfy stdDev > 0." << nl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::randomVariable::randomVariable()
:
    type_(),
    weights_(),
    values_(),
    lowerBound_(-GREAT),
    upperBound_(GREAT),
    meanValue_(),
    standardDev_(),
    maxValue_(),
    minValue_(),
    mu_(),
    sigma_(),
    generateCounter_(0)
{
    initializeRandomGenerator();
}


Foam::randomVariable::randomVariable
(
    const dictionary& dict
)
:
    type_(),
    weights_(),
    values_(),
    lowerBound_(-GREAT),
    upperBound_(GREAT),
    meanValue_(),
    standardDev_(),
    maxValue_(),
    minValue_(),
    mu_(),
    sigma_(),
    generateCounter_(0)
{
    read(dict);
    initializeRandomGenerator();
}


Foam::randomVariable::randomVariable
(
    const dictionary& dict,
    const scalar lowerBound,
    const scalar upperBound
)
:
    type_(),
    weights_(),
    values_(),
    lowerBound_(-lowerBound),
    upperBound_(upperBound),
    meanValue_(),
    standardDev_(),
    maxValue_(),
    minValue_(),
    mu_(),
    sigma_(),
    generateCounter_(0)
{
    read(dict);
    initializeRandomGenerator();
}


Foam::randomVariable::randomVariable
(
    Istream& is
)
:
    type_(),
    weights_(),
    values_(),
    lowerBound_(-GREAT),
    upperBound_(GREAT),
    meanValue_(),
    standardDev_(),
    maxValue_(),
    minValue_(),
    mu_(),
    sigma_(),
    generateCounter_(0)
{
    dictionary dict(is);
    read(dict);
    initializeRandomGenerator();
}


Foam::randomVariable::randomVariable
(
    const randomVariable& obj
)
:
    type_(obj.type_),
    weights_(obj.weights_),
    values_(obj.values_),
    lowerBound_(obj.lowerBound_),
    upperBound_(obj.upperBound_),
    meanValue_(obj.meanValue_),
    standardDev_(obj.standardDev_),
    maxValue_(obj.maxValue_),
    minValue_(obj.minValue_),
    mu_(obj.mu_),
    sigma_(obj.sigma_),
    generateCounter_(0)
{
    initializeRandomGenerator();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::randomVariable::~randomVariable()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar
Foam::randomVariable::meanValue() const
{
    return meanValue_;
}


scalar
Foam::randomVariable::generate() const
{
    scalar val;
    scalar spacing;
    generateCounter_++;

    switch (type_)
    {
        case CONSTANT:
            return meanValue_;

        case CYCLIC_LIST:
            return values_[(generateCounter_ - 1) % values_.size()];

        case DISCRETE:
            int i;
            i = 0;
            scalar rightBinSide;
            rightBinSide = weights_[i];
            val = scalar01();
            while (val > rightBinSide)
            {
                i++;
                rightBinSide += weights_[i];
            }
            return values_[i];
        
        case UNIFORM:
            return minValue_ + scalar01()*(maxValue_ - minValue_);

        case GAUSSIAN:
            val = meanValue_ + GaussNormal()*standardDev_;
            return max(min(maxValue_, val), minValue_);

        case LOG_NORMAL_SPACING:
            // generate a log-normal relative spacing
            spacing = exp(mu_ + GaussNormal()*sigma_);
            // transform it to a line density
            return 1./(1. + spacing);
    }

    return -1;
}


void randomVariable::write(Ostream& os) const
{
    os.writeKeyword("type")        << type_        << token::END_STATEMENT;
    os.writeKeyword("lowerBound")  << lowerBound_  << token::END_STATEMENT;
    os.writeKeyword("upperBound")  << upperBound_  << token::END_STATEMENT;
    os.writeKeyword("meanValue")   << meanValue_   << token::END_STATEMENT;
    os.writeKeyword("standardDev") << standardDev_ << token::END_STATEMENT;
    os.writeKeyword("maxValue")    << maxValue_    << token::END_STATEMENT;
    os.writeKeyword("minValue")    << minValue_    << token::END_STATEMENT;
    os.writeKeyword("mu")          << mu_          << token::END_STATEMENT;
    os.writeKeyword("sigma")       << sigma_       << token::END_STATEMENT;
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Istream&
Foam::operator>>(Istream& is, randomVariable& obj)
{
    dictionary dict(is);
    obj.read(dict);

    return is;
}


Foam::Ostream&
Foam::operator<<(Ostream& os, const randomVariable& obj)
{
    obj.write(os);

    return os;
}

// ************************************************************************* //
