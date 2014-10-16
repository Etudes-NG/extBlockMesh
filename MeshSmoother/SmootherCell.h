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

#ifndef MESHSMOOTHERCELL_H
#define MESHSMOOTHERCELL_H

#include "cellShape.H"
#include "polyMesh.H"

#include "SmootherBoundary.h"
#include "SmootherPoint.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                      Class MeshSmootherCell Declaration
\*---------------------------------------------------------------------------*/

class SmootherCell
{

    //- Private data

        // Pointer of smoother boundary
        static SmootherBoundary* _bnd;

        // Transformation parameter
        static scalar _transParam;

        // Cell quality (mean ratio)
        scalar _quality;

        // Reference of cell shape
        const cellShape& _cellShape;

    //- Private member functions

        // Tetrahedral quality computation
        scalar tetCellQuality(const label ref) const;

        // get point from cell position
        inline const point& initPt(const label p) const;
        inline const point& relaxPt(const label p) const;

        scalar fastPow(const scalar& s) const;

public:

    //- Constructors

        //- Construct from cellShape
        SmootherCell(const cellShape& cell);


    //- Member functions

        // Set/get cell quality
        const scalar& quality() const{return _quality;}
        void computeQuality();

        // Transform cell
        pointField geometricTransform();

        // Set/get polyMesh
        void setStaticItems(SmootherBoundary *bnd, const scalar& t);
};

const point &SmootherCell::initPt(const label p) const
{
    return _bnd->pt(_cellShape[p])->getInitialPoint();
}

const point &SmootherCell::relaxPt(const label p) const
{
    return _bnd->pt(_cellShape[p])->getRelaxedPoint();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // MESHSMOOTHERCELL_H

// ************************************************************************* //
