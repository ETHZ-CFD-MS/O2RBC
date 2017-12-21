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
    Foam::polygonalTube

Description
    Circular tube around a polygonal edge. The diameters are considered to be
    constant along the each edge segment.

SourceFiles
    polygonalTube.C

\*---------------------------------------------------------------------------*/

#ifndef polygonalTube_H
#define polygonalTube_H

#include "circularTube.H"

#include "polygonalEdge.H"
#include "ellipseAxes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;

// Forward declaration of friend functions and operators
class polygonalTube;
Ostream& operator<<(Ostream&, const polygonalTube&);


/*---------------------------------------------------------------------------*\
                         Class polygonalTube Declaration
\*---------------------------------------------------------------------------*/

class polygonalTube
: 
    public circularTube
{
    // Private data
    
        //- Reference to a polygonalEdge.
        const polygonalEdge& edge_;

        //- List of edge segment inner diameters
        scalarList innerDiameters_;

        //- Ratio between outer and inner tube diameter
        scalar outerInnerDiameterRatio_;

        //- List of ellipseAxes at vertices of the polygonal edge path.
        List<ellipseAxes> axesList_;

    // Private Member Functions

        //- Check and adapt diameters
        void checkAndAdaptDiameters
        (
            const dictionary& tubeOptions
        );

        //- Compute the list of ellispse axes.
        void computeEllipseAxesList();

        //- Convert the segment diameters to the effective diameter.
        //  A call to this method overwrites the segment diameters!
        void convertToEffectiveDiameter();

        //- Disallow default bitwise copy construct
        polygonalTube(const polygonalTube&);

        //- Disallow default bitwise assignment
        void operator=(const polygonalTube&);

public:

    // Constructors

        //- Construct from components
        polygonalTube
        (
            const polygonalEdge& edge,
            const scalarList& innerDiameters,
            const ellipseAxes& startAxes,
            const ellipseAxes& endAxes,
            const dictionary& tubeOptions
        );

        //- Construct from components, without specifying ellipseAxes
        polygonalTube
        (
            const polygonalEdge& edge,
            const scalarList& innerDiameters,
            const dictionary& tubeOptions
        );

        //- Construct with a single diameter
        polygonalTube
        (
            const polygonalEdge& edge,
            const scalar innerDiameter,
            const ellipseAxes& startAxes,
            const ellipseAxes& endAxes,
            const scalar outerInnerDiameterRatio
        );

    //- Destructor
    virtual ~polygonalTube();


    // Member Functions

        // Access

            const polygonalEdge& edge() const
            {
                return edge_;
            }

            const scalarList segmentOuterDiameters() const;
            const scalarList segmentInnerDiameters() const;

        // Setters
            
            virtual void setStartAxes(const ellipseAxes&);
            virtual void setEndAxes(const ellipseAxes&);

        // Geometric queries

            //- Mean outer tube diameter
            virtual scalar meanOuterDiameter() const;

            //- Mean inner tube diameter
            virtual scalar meanInnerDiameter() const;

            //- Return the tube outer diameter at a given coordinate.
            virtual scalar outerDiameter(const scalar s) const;
    
            //- Return the tube inner diameter at a given coordinate.
            virtual scalar innerDiameter(const scalar s) const;
    
            //- Compute the effective diameter of the tube.
            //  The effective diameter yields the same tube resistance as the
            //  sequence of segment diameters.
            scalar effectiveDiameter() const;

            //- Return the closest ellipseAxes with coordinate < s
            virtual ellipseAxes ellipseAxesBefore(const scalar s) const;

            //- Return the closest ellipseAxes with coordinate >= s
            virtual ellipseAxes ellipseAxesAfter(const scalar s) const;

            //- Return the ellipse on the start side of the given segment.
            ellipseAxes ellipseAxesAtSegmentStart(const label segmentI) const;

            //- Return the ellipse on the end side of the given segment.
            ellipseAxes ellipseAxesAtSegmentEnd(const label segmentI) const;
            
            //- Return a skeleton of the point tubes
            virtual pointField skeletonPoints() const;

        // Visitor call
            
            //- Call the visitor's visit method (see Visitor pattern).
            virtual void accept(tubeVisitor&) const;

        // Write

            virtual void write(Ostream& os) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const polygonalTube&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
