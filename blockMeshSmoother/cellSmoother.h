/*---------------------------------------------------------------------------*\
  extBlockMesh
  Copyright (C) 2014 Etudes-NG
  ---------------------------------
License
    This file is part of extBlockMesh.

    extBlockMesh is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    extBlockMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with extBlockMesh.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef CELLSMOOTHER_H
#define CELLSMOOTHER_H

#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cellSmoother Declaration
\*---------------------------------------------------------------------------*/

class cellSmoother
{
    // Private data

        // Cell points
        pointField &points_;

    // Private member functions

        //- tetrahedral mean ratio
        scalar tetrahedralMeanRatio
        (
            const label &pt,
            const label &pt1,
            const label &pt2,
            const label &pt3
        ) const;

        //- Hexahedron dual octahedron
        pointField dualOctahedron() const;

        //- Octahedron faces centroid
        pointField dualOctahedronFaceCentroid(const pointField &oct) const;

        //- Octahedron normal
        pointField dualOctahedronNormals(const pointField &oct) const;

        //- Transfomred hex
        pointField tranformedHexahedron
        (
            const scalar &cor,
            const pointField &octC,
            const pointField &octN
        ) const;

        //- Centroid of hex
        point centroidOfHex(const pointField &Hp) const;

        //- Ratio of length
        scalar ratioOfAvgLength(pointField &Hp) const;

        //- Edge average length
        scalar edgeAverageLength() const;

public:
    // Constructors

        //- Construct from blockMesh and dictionary
        cellSmoother(pointField &H);

    // Member functions

        //- meanRatio
        scalar meanRatio() const;

        // geometricTranform
        // FIXME decompose geometricTransfom in multiple methods (trans, scal,.)
        pointField geometricTranform(const scalar &cor) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // CELLSMOOTHER_H

// ************************************************************************* //
