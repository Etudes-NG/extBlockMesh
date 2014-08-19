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

#ifndef POINTTOPO_H
#define POINTTOPO_H

#include "blockMesh.H"

#include <map>
#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockMeshTopology;

/*---------------------------------------------------------------------------*\
             Class pointTopo Declaration
\*---------------------------------------------------------------------------*/

class pointTopo
{
    //- Private data

    //- Private member functions

protected:
        //- Pointer of boundary points topology
        blockMeshTopology *topo_;

    //- Protected member functions


public:
    //- Constructors

        //- Construct from
        pointTopo(blockMeshTopology *topo);

    //- Destructor
        virtual ~pointTopo() = 0;

    //- Member functions

        //- Point smoothed with respect of topology constraints
        virtual point smoothedPoint(const point &guessedPoint);

        virtual std::map<scalar,point> mapNeiborFeaturePts
        (const point &guessedPoint) const;

        //- Get linked point of feature edge point
        virtual std::set<label> getPointLinked() const;

        //- Get linked triangles faces of feature edge point
        virtual std::set<std::set<label> > getTrianglesLinked() const;


        //- Get optimal point from guessed
        virtual point getFeatureEdgePoint
        (
            const point &guessedPoint,
            const label &ref
        );

        //- Get boundaryPoint coord
        virtual point &getboundaryPoint();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // POINTTOPO_H

// ************************************************************************* //


