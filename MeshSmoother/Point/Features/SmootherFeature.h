#ifndef MESHSMOOTHERPOINTFEATURE_H
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

#define MESHSMOOTHERPOINTFEATURE_H

#include "SmootherPoint.h"

#include "SmootherBoundary.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class MeshSmootherPointSurface Declaration
\*---------------------------------------------------------------------------*/

class SmootherFeature
:
    public SmootherPoint
{
    //- Private data

protected:

    //- Protected data

        // Feature ref
        label _featureRef;

public:
    SmootherFeature(const label ref, const label featureRef);

    ~SmootherFeature() {}

    inline void needSnap();
};

void SmootherFeature::needSnap()
{
    if (_relaxLevel == 0)
    {
        _bnd->removeSnapPoint(_ptRef);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // MESHSMOOTHERPOINTFEATURE_H

// ************************************************************************* //
