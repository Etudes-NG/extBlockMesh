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

#ifndef FEATUREEDGEPOINT_H
#define FEATUREEDGEPOINT_H

#include "boundaryPoint.h"

#include <set>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class pointTopo Declaration
\*---------------------------------------------------------------------------*/

class featureEdgePoint : public boundaryPoint
{
    //- Private data

        //- Point link at beginning
        std::set<label> pointLinked_;

        //- Point link during move
        std::set<label> pntLinkedNew_;

        //- Initial point label
        label featurePtLabel_;

        //- Point label during smoothing
        label newFeaturePtLabel_;

        // TODO change "new" by pointer of pointTopo

    //- Private member functions

        //- Get linked point of feature edge point
        std::set<label> getPointLinked() const;

public:
    //- Constructors

        //- Construct from
        featureEdgePoint
        (
            const std::set<label> &pointLinked,
            const std::set<std::set<Foam::label> > &triangles,
            const point &initialPoint,
            const label &initialLabel,
            blockMeshTopology *topo
        );

    //- Destructor
        ~featureEdgePoint();

    //- Member functions

        //- Get optimal point with repect of point topo
        point smoothedPoint(const point &guessedPoint);

        //- Get map of points sorted by minimal distance
        std::map<scalar,point> mapNeiborFeaturePts
        (
            const point &guessedPoint
        ) const;

        //- Change feature edge linked points and get optimal point
        point changeFeatureEdgeLinkedsPoint
        (
            const label &newRef,
            const point &guessedPoint
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // FEATUREEDGEPOINT_H

// ************************************************************************* //

