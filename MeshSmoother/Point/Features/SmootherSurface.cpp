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

#include "SmootherSurface.h"

#include "polyMesh.H"

#include "SmootherBoundary.h"

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SmootherSurface::SmootherSurface(const label ref, const label featureRef)
:
    SmootherFeature(ref, featureRef)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SmootherSurface::GETMeSmooth()
{
    SmootherPoint::GETMeSmooth();
    _movedPt = _bnd->snapToSurf(_featureRef, _movedPt);
}

void SmootherSurface::snap()
{
    _movedPt = _bnd->snapToSurf(_featureRef, _initialPt);
}

void SmootherSurface::featLaplaceSmooth()
{
    const labelList& pp = _polyMesh->pointPoints(_ptRef);
    label nbPt = 0;
    _movedPt = point(0.0, 0.0, 0.0);
    forAll(pp, ptI)
    {
        if (_bnd->pt(pp[ptI])->isSurface())
        {
            _movedPt += _bnd->pt(pp[ptI])->getRelaxedPoint();
            ++nbPt;
        }
    }
    _movedPt /= nbPt;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
