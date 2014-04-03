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

#include "featureEdgePoint.h"

#include "blockMeshTopology.h"

#include <algorithm>
#include <map>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::featureEdgePoint::featureEdgePoint
(
    const std::set<label> &pointLinked,
    const std::set<std::set<Foam::label> > &triangles,
    const point &initialPoint,
    const label &initialLabel,
    blockMeshTopology *topo
)
    :
    boundaryPoint(triangles, initialPoint, initialLabel, topo),
    pointLinked_(pointLinked),
    pntLinkedNew_(pointLinked),
    featurePtLabel_(initialLabel),
    newFeaturePtLabel_(initialLabel)
{
}

Foam::featureEdgePoint::~featureEdgePoint()
{

}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

std::map<Foam::scalar, Foam::point> Foam::featureEdgePoint::mapNeiborFeaturePts
(
    const point &guessedPoint
) const
{
    std::map<scalar,point> minDist;

    for
    (
        std::set<label>::iterator ptI = pntLinkedNew_.begin();
        ptI != pntLinkedNew_.end();
        ++ptI
    )
    {
        const point &p1(topo_->getBndPointCoord(*ptI));
        const point &p2(topo_->getBndPointCoord(newFeaturePtLabel_));
        const point &p3(guessedPoint);
        const scalar u
        (
            (
                (p3.x() - p1.x())*(p2.x() - p1.x()) +
                (p3.y() - p1.y())*(p2.y() - p1.y()) +
                (p3.z() - p1.z())*(p2.z() - p1.z())
            ) /
            (magSqr(p2 - p1) + VSMALL)
        );
        const point px
        (
            p1.x() + u*(p2.x() - p1.x()),
            p1.y() + u*(p2.y() - p1.y()),
            p1.z() + u*(p2.z() - p1.z())
        );
        if (u > 0.0 && u <= 1.0)
        { // closest point fall within the line segment
            minDist.insert(std::make_pair<scalar,point>(mag(p3 - px), px));
        }
    }

    return minDist;
}

std::set<Foam::label> Foam::featureEdgePoint::getPointLinked() const
{
    return pointLinked_;
}

Foam::point Foam::featureEdgePoint::changeFeatureEdgeLinkedsPoint
(
    const label &newRef,
    const point &guessedPoint
)
{
    // Update the points linked
    pntLinkedNew_ = topo_->getPointTopoPtr(newRef)->getPointLinked();
    newFeaturePtLabel_ = newRef;

    return smoothedPoint(guessedPoint);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::featureEdgePoint::smoothedPoint(const point &guessedPoint)
{
    std::map<scalar,point> minDists(mapNeiborFeaturePts(guessedPoint));

    if (minDists.empty())
    {
        const scalar distCenter
        (
            mag(guessedPoint - topo_->getBndPointCoord(newFeaturePtLabel_))
        );
        const scalar distExtrem1
        (
            mag(guessedPoint - topo_->getBndPointCoord(*pntLinkedNew_.begin()))
        );
        const scalar distExtrem2
        (
            mag(guessedPoint - topo_->getBndPointCoord(*pntLinkedNew_.rbegin()))
        );

        if (distCenter < distExtrem1 && distCenter < distExtrem2)
        { // Boundary is convex
            return topo_->getBndPointCoord(newFeaturePtLabel_);
        }
        else if (distExtrem1 < distExtrem2)
        { // Nearest boundary point is begin

            return changeFeatureEdgeLinkedsPoint
            (
                *pntLinkedNew_.begin(),
                guessedPoint
            );
        }
        else
        { // Nearest boundary point is rbegin

            return changeFeatureEdgeLinkedsPoint
            (
                *pntLinkedNew_.rbegin(),
                guessedPoint
            );
        }
    }
    else
    {
        return minDists.begin()->second;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
