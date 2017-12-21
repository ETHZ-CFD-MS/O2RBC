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
    Foam::piecewiseConstantInterpolationTable

Description
    An interpolation/look-up table of scalar vs \<Type\> values.
    The reference scalar values must be monotonically increasing.

    The interpolation between values is piecewise constant. For t1 <= t < t2,
    the interpolated value for t is the value at time t1.

    The handling of out-of-bounds values depends on the current setting
    of \a outOfBounds.

    If \a REPEAT is chosen for the out-of-bounds handling, the final time
    value is treated as being equivalent to time=0 for the following periods.

    The code was taken from interpolationTable. Inheritance could not be
    used since interpolationTable does not support Type=List<T>.


Note
    - Accessing an empty list results in an error.
    - Accessing a list with a single element always returns the same value.

SourceFiles
    piecewiseConstantInterpolationTable.C

\*---------------------------------------------------------------------------*/

#ifndef piecewiseConstantInterpolationTable_H
#define piecewiseConstantInterpolationTable_H

#include "List.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class piecewiseConstantInterpolationTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class piecewiseConstantInterpolationTable
:
    public List<Tuple2<scalar, Type> >
{
public:

    // Public data types

        //- Enumeration for handling out-of-bound values
        enum boundsHandling
        {
            ERROR,          /*!< Exit with a FatalError */
            WARN,           /*!< Issue warning and clamp value (default) */
            CLAMP,          /*!< Clamp value to the start/end value */
            REPEAT          /*!< Treat as a repeating list */
        };


private:

    // Private data

        //- Enumeration for handling out-of-bound values
        boundsHandling boundsHandling_;

        //- File name
        fileName fileName_;

    // Private Member Functions

        //- Read the table of data from file
        void readTable();


public:

    // Constructors

        //- Construct null
        piecewiseConstantInterpolationTable();

        //- Construct from components
        piecewiseConstantInterpolationTable
        (
            const List<Tuple2<scalar, Type> >& values,
            const typename piecewiseConstantInterpolationTable<Type>::boundsHandling bounds,
            const fileName& fName
        );

        //- Construct given the name of the file containing the table of data
        piecewiseConstantInterpolationTable(const fileName& fName);

        //- Construct by reading the data from a dictionary with a given key.
        piecewiseConstantInterpolationTable(const dictionary& dict, const word& key);

        //- Construct by reading the fileName and boundsHandling from dictionary
        //  and read the table from that file.
        //  This is a specialised constructor used by patchFields
        piecewiseConstantInterpolationTable(const dictionary& dict);

        //- Construct copy
        piecewiseConstantInterpolationTable
        (
            const piecewiseConstantInterpolationTable& interpTable
        );


    // Member Functions

        //- Return the out-of-bounds handling as a word
        word boundsHandlingToWord(const boundsHandling& bound) const;

        //- Return the out-of-bounds handling as an enumeration
        boundsHandling wordToBoundsHandling(const word& bound) const;

        //- Set the out-of-bounds handling from enum, return previous setting
        boundsHandling outOfBounds(const boundsHandling& bound);

        //- Check that list is monotonically increasing
        //  Exit with a FatalError if there is a problem
        void check() const;

        //- Write
        void write(Ostream& os) const;

        //- Return the rate of change at the interpolation location
        //  for the given value
        Type rateOfChange(const scalar) const;

        //- Return whether the value changed between the two given interpolation
        //  locations.
        bool valueChanged(const scalar, const scalar) const;


    // Member Operators

        //- Return an element of constant Tuple2<scalar, Type>
        const Tuple2<scalar, Type>& operator[](const label) const;

        //- Return a value interpolated using a piecewise constant approximation
        Type operator()(const scalar) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "piecewiseConstantInterpolationTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
