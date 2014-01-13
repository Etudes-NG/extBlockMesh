#include "featureEdgePoint.h"

#include "blockMeshTopology.h"

#include <algorithm>
#include <map>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::featureEdgePoint::featureEdgePoint
(
    const std::set<label> &pointLinked,
    blockMeshTopology *topo
)
    :
      pointLinked_(pointLinked),
      pointLinkedNew_(pointLinked),
      topo_(topo)
{
}

Foam::featureEdgePoint::~featureEdgePoint()
{

}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

std::map<Foam::scalar, Foam::point> Foam::featureEdgePoint::minDist
(
    const point &guessedPoint,
    const label &pointRef
)
{
    std::map<scalar,point> minDist;

    for
    (
        std::set<label>::iterator ptI = pointLinkedNew_.begin();
        ptI != pointLinkedNew_.end();
        ++ptI
    )
    {
        const point p1(topo_->getBndPt(*ptI));
        const point p2(topo_->getBndPt(pointRef));
        const point p3(guessedPoint);
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

Foam::point Foam::featureEdgePoint::getPoint
(
    const point &guessedPoint,
    const label &ref
)
{
    std::map<scalar,point> minDists(minDist(guessedPoint, ref));

    if (minDists.empty())
    {
        const scalar distCenter
        (
            mag(guessedPoint - topo_->getBndPt(ref))
        );
        const scalar distExtrem1
        (
            mag(guessedPoint - topo_->getBndPt(*pointLinkedNew_.begin()))
        );
        const scalar distExtrem2
        (
            mag(guessedPoint - topo_->getBndPt(*pointLinkedNew_.rbegin()))
        );

        if (distCenter < distExtrem1 && distCenter < distExtrem2)
        { // Boundary is convex
            return topo_->getBndPt(ref);
        }
        else if (distExtrem1 < distExtrem2)
        { // Nearest boundary point is begin

            return changeLinkedsPoint(*pointLinkedNew_.begin(), guessedPoint);
        }
        else
        { // Nearest boundary point is rbegin

            return changeLinkedsPoint(*pointLinkedNew_.rbegin(), guessedPoint);
        }
    }
    else
    {
        return minDists.begin()->second;
    }
}

Foam::point Foam::featureEdgePoint::changeLinkedsPoint
(
    const label &newRef,
    const point &guessedPoint
)
{
    pointLinkedNew_ = topo_->getPointTopoPtr(newRef)->getPointLinked();

    if (pointLinkedNew_.size() == 2)
    { // search point with this new linked point
        return getPoint(guessedPoint, newRef);
    }
    else
    { // Point linked is extremity
        Info<< "Error feature edge point next invalid\n";
        return topo_->getBndPt(newRef);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::featureEdgePoint::smoothedPoint
(
    const point &guessedPoint,
    const label &pointRef
)
{
    if (pointLinked_.size() == 2)
    { // Point is linked to 2 points
        return getPoint(guessedPoint, pointRef);
    }
    else if (pointLinked_.size() > 2 || pointLinked_.size() == 1)
    { // Corner point
        return topo_->getBndPt(pointRef);
        Info<< "Corner point\n";
    }
    else
    {
        Info<< "Error feature edge point\n";
        return topo_->getBndPt(pointRef);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
