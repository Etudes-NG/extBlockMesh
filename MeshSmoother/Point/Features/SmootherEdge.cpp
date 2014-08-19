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

#include "SmootherEdge.h"

#include "SmootherBoundary.h"

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SmootherEdge::SmootherEdge(const label ref, const label featureRef)
:
    SmootherFeature(ref, featureRef)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SmootherEdge::GETMeSmooth()
{
    SmootherPoint::GETMeSmooth();
    _movedPt = _bnd->snapToEdge(_featureRef, _movedPt);
}

void SmootherEdge::snap()
{
    const labelList& pp = _polyMesh->pointPoints(_ptRef);
    label nbPt = 0;
    _movedPt = point(0.0, 0.0, 0.0);
    forAll(pp, ptI)
    {
        if (_bnd->pt(pp[ptI])->isEdge())
        {
            _movedPt += _bnd->pt(pp[ptI])->getInitialPoint();
            ++nbPt;
        }
    }
    _movedPt /= nbPt;

    _movedPt = _bnd->snapToEdge(_featureRef, _movedPt);
}

void SmootherEdge::featLaplaceSmooth()
{
    const labelList& pp = _polyMesh->pointPoints(_ptRef);
    label nbPt = 0;
    _movedPt = point(0.0, 0.0, 0.0);
    forAll(pp, ptI)
    {
        if (_bnd->pt(pp[ptI])->isEdge())
        {
            _movedPt += _bnd->pt(pp[ptI])->getRelaxedPoint();
            ++nbPt;
        }
    }

    _movedPt /= nbPt;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
