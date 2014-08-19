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

#include "SmootherPoint.h"

#include "polyMesh.H"

#include "SmootherParameter.h"
#include "SmootherBoundary.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    polyMesh* SmootherPoint::_polyMesh = NULL;
    SmootherBoundary* SmootherPoint::_bnd = NULL;
    SmootherParameter* SmootherPoint::_param = NULL;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SmootherPoint::SmootherPoint(const label ref)
:
    _relaxedPt(_polyMesh->points()[ref]),
    _ptRef(ref)
{
}

Foam::SmootherPoint::SmootherPoint()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SmootherPoint::setStaticItems
(
    SmootherBoundary *bnd,
    SmootherParameter *param,
    polyMesh *mesh
)
{
    _bnd = bnd;
    _param = param;
    _polyMesh = mesh;
}



void SmootherPoint::laplaceSmooth()
{
    const labelList& pp = _polyMesh->pointPoints(_ptRef);
    _movedPt = point(0.0, 0.0, 0.0);
    forAll(pp, ptI)
    {
        _movedPt += _bnd->pt(pp[ptI])->getInitialPoint();
    }
    _movedPt /= pp.size();
}





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
