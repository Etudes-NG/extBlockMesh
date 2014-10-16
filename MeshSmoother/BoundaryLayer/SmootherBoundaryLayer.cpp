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

#include "SmootherBoundaryLayer.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SmootherBoundaryLayer::SmootherBoundaryLayer(const dictionary &blDict)
:
    _nbLayers(readLabel(blDict.lookup("nSurfaceLayers"))),
    _expansionRatio(readScalar(blDict.lookup("expansionRatio"))),
    _relativeSize(readBool(blDict.lookup("relativeSizes"))),
    _finalLayerThickness(readScalar(blDict.lookup("finalLayerThickness")))
{
    Info<< "        - Number of BL       : " << _nbLayers << nl;
    Info<< "        - Expansion ratio    : " << _expansionRatio << nl;
    Info<< "        - Relatice size      : " << _relativeSize << nl;
    Info<< "        - Final thickness    : " << _finalLayerThickness << nl;
}

Foam::SmootherBoundaryLayer::SmootherBoundaryLayer()
:
    _nbLayers(0.0),
    _expansionRatio(1.0),
    _relativeSize(true),
    _finalLayerThickness(1)
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
