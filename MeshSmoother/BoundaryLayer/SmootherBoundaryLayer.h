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

#ifndef SMOOTHERBOUNDARYLAYER_H
#define SMOOTHERBOUNDARYLAYER_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class SmootherBoundaryLayer Declaration
\*---------------------------------------------------------------------------*/

class SmootherBoundaryLayer
{
    //- Private data

        // Number of boundary layers
        label _nbLayers;

        // Expansion ratio
        scalar _expansionRatio;

        // Relative size
        bool _relativeSize;

        // Wanted thickness of final added cell layer. If multiple
        // layers is the thickness of the layer furthest away from the
        // wall. Relative to undistorted size of cell outside layer.
        scalar _finalLayerThickness;

public:
    SmootherBoundaryLayer(const dictionary &blDict);
    SmootherBoundaryLayer();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // SMOOTHERBOUNDARYLAYER_H

// ************************************************************************* //
