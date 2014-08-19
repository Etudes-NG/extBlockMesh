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

#include "cornerPoint.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cornerPoint::cornerPoint
(
    const std::set<std::set<Foam::label> > &triangles,
    const point &initialPoint,
    const label &initialLabel,
    blockMeshTopology *topo
)
    :
      boundaryPoint(triangles, initialPoint, initialLabel, topo)
{
}

Foam::cornerPoint::~cornerPoint()
{

}


// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::cornerPoint::smoothedPoint
(
    const point &guessedPoint,
    const label &pointRef
) const
{
    // No change as the point is fixed
    return initialPoint_;
}

Foam::point Foam::cornerPoint::changeFeatureEdgeLinkedsPoint
(
    const Foam::label &newRef,
    const Foam::point &guessedPoint
)
{
    FatalErrorIn("changeFeatureEdgeLinkedsPoint(guessedPoint, ref)")
        << "Accessed from a non feature edge point\n"
        << nl
        << exit(FatalError);

    return point();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //

