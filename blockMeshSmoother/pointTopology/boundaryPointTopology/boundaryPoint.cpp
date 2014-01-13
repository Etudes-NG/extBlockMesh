#include "boundaryPoint.h"

#include "blockMeshTopology.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryPoint::boundaryPoint
(
    const std::set<std::set<Foam::label> > &triangles,
    blockMeshTopology *topo
)
    :
      pointTopo(triangles, topo)
{
}

Foam::boundaryPoint::~boundaryPoint()
{
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

bool Foam::boundaryPoint::isOnTriangle
(
    const label &refP2,
    const label &refP3,
    const point &p,
    point &out,
    const label &ref
) const
{
    // http://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle/

    const point p1(topo_->getBndPt(ref));
    const point p2(topo_->getBndPt(refP2));
    const point p3(topo_->getBndPt(refP3));

    const point u(p2 - p1);
    const point v(p3 - p1);
    const point n(u ^ v);
    const point w(p - p1);

    const scalar lambda(fabs((u ^ w) & n)/(n & n));
    const scalar beta(fabs((w ^ v) & n)/(n & n));
    const scalar alpha(1.0 - lambda - beta);

    if
    (
        alpha >= 0 && alpha <= 1 &&
        lambda >= 0 && lambda <= 1 &&
        beta >= 0 && beta <= 1
    )
    { // Point is on triangle
        out = alpha*p1 + beta*p2 + lambda*p3;
        return true;
    }
    else
    {
        return false;
    }
}

std::set<std::set<Foam::label> > Foam::boundaryPoint::getTrianglesLinked() const
{
    return triangles_;
}

Foam::point Foam::boundaryPoint::projectedBndPoint
(
    const point &guessedPoint,
    const label &ref
)
{
    std::map<scalar,point> minDists(minDist(guessedPoint, ref));

    if (minDists.empty())
    {
        const scalar distCenter(mag(guessedPoint - topo_->getBndPt(ref)));

        // Set of all extremity point
        std::set<label> extremPoint;
        for // all triangle
        (
            std::set<std::set<label> >::iterator triI = trianglesNew_.begin();
            triI != trianglesNew_.end();
            ++triI
        )
        {
            for // all point
            (
                std::set<label>::iterator ptI = triI->begin();
                ptI != triI->end();
                ++ptI
            )
            {
                 extremPoint.insert(*ptI);
            }
        }

        // Store the dist between guessed point and extremity point
        std::map<scalar, label> dist;

        for // all extremity point
        (
            std::set<label>::iterator ptI = extremPoint.begin();
            ptI != extremPoint.end();
            ++ptI
        )
        {
            dist.insert
            (
                std::make_pair<scalar,point>
                (
                    mag(guessedPoint - topo_->getBndPt(*ptI)),
                    *ptI
                )
            );
        }

        if (distCenter < dist.begin()->first)
        { // convex
            return topo_->getBndPt(ref);
        }
        else
        {
            return changeLinkedsPoint(dist.begin()->second, guessedPoint);
        }
    }
    else
    {
        return minDists.begin()->second;
    }
}

Foam::point Foam::boundaryPoint::changeLinkedsPoint
(
    const Foam::label &newRef,
    const Foam::point &guessedPoint
)
{
    return topo_->getPointTopoPtr(newRef)->getPoint(guessedPoint, newRef);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::point Foam::boundaryPoint::smoothedPoint
(
    const Foam::point &guessedPoint,
    const label &pointRef
)
{
    point pt;
    for
    (
        std::set<std::set<label> >::iterator triI = trianglesNew_.begin();
        triI != trianglesNew_.end();
        ++triI
    )
    {
        if
        (
             isOnTriangle
             (
                 *(*triI).begin(),
                 *(*triI).rbegin(),
                 guessedPoint,
                 pt,
                 pointRef
             )
         )
        {
            return pt;
        }
    }



    Info<< "Error no new triangle found\n";

    // TODO: compute guessed point from normal and topo
    return guessedPoint;
}

std::map<scalar, Foam::point> Foam::boundaryPoint::minDist
(
    const point &guessedPoint,
    const label &pointRef
)
{
    std::map<scalar,point> minDist;

    for
    (
        std::set<std::set<label> >::iterator triI = trianglesNew_.begin();
        triI != trianglesNew_.end();
        ++triI
    )
    {
        point pt;
        if
        (
             isOnTriangle
             (
                 *(*triI).begin(),
                 *(*triI).rbegin(),
                 guessedPoint,
                 pt,
                 pointRef
             )
         )
        {
            minDist.insert
            (
                std::make_pair<scalar,point>(mag(guessedPoint - pt), pt)
            );
        }
    }

    return minDist;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
