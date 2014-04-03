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

#include "pointTopo.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointTopo::pointTopo(blockMeshTopology *topo)
    :
      topo_(topo)
{
}

Foam::pointTopo::~pointTopo()
{

}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::set<Foam::label> Foam::pointTopo::getPointLinked() const
{
    FatalErrorIn("getPointLinked()")
        << "Accessed from a non feature edge point\n"
        << nl
        << exit(FatalError);

    return std::set<label>();
}

Foam::point Foam::pointTopo::smoothedPoint(const Foam::point &guessedPoint)
{
    FatalErrorIn("smoothedPoint()")
        << "Accessed from a non feature edge point\n"
        << nl
        << exit(FatalError);

    return point();
}

std::map<Foam::scalar, Foam::point> Foam::pointTopo::mapNeiborFeaturePts
(
    const Foam::point &guessedPoint
) const
{
    FatalErrorIn("mapNeiborFeaturePts()")
        << "Accessed from a non feature edge point\n"
        << nl
        << exit(FatalError);


    return std::map<scalar,point>();
}

std::set<std::set<Foam::label> > Foam::pointTopo::getTrianglesLinked() const
{
    FatalErrorIn("getTrianglesLinked()")
        << "Accessed from a non feature edge point\n"
        << nl
        << exit(FatalError);

    return std::set<std::set<Foam::label> >();
}

Foam::point Foam::pointTopo::getFeatureEdgePoint
(
    const Foam::point &guessedPoint,
    const Foam::label &ref
)
{
    FatalErrorIn("getFeatureEdgePoint()")
        << "Accessed from a non feature edge point\n"
        << nl
        << exit(FatalError);

    return point();
}

Foam::point &Foam::pointTopo::getboundaryPoint()
{
    FatalErrorIn("getboundaryPoint()")
        << "Accessed from a non boundary point\n"
        << nl
        << exit(FatalError);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
