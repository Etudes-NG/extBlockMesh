#include "featureEdgePoint.h"

#include "blockMeshTopology.h"

#include <algorithm>
#include <map>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::featureEdgePoint::featureEdgePoint()
{
}

Foam::featureEdgePoint::~featureEdgePoint()
{

}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::featureEdgePoint::smoothedPoint
(
    const point &guessedPoint,
    const blockMesh *blocks,
    const label &pointRef,
    const blockMeshTopology *topo
) const
{
    if (pointLinked_.size() > 2)
    { // Corner point
        return topo->getBoundaryPoint(pointRef);
        Info<< "Corner point\n";
    }
    else if (pointLinked_.size() == 2)
    {
        std::map<scalar,point> minDist;

        for
        (
            std::set<label>::iterator ptI = pointLinked_.begin();
            ptI != pointLinked_.end();
            ++ptI
        )
        {
            const point p1(topo->getBoundaryPoint(*ptI));
            const point p2(topo->getBoundaryPoint(pointRef));
            const point p3(guessedPoint);
            const scalar u
            (
                (
                    (p3.x() - p1.x())*(p2.x() - p1.x()) +
                    (p3.y() - p1.y())*(p2.y() - p1.y()) +
                    (p3.z() - p1.z())*(p2.z() - p1.z())
                ) /
                magSqr(p2 - p1)
            );
            const point px
            (
                p1.x() + u*(p2.x() - p1.x()),
                p1.y() + u*(p2.y() - p1.y()),
                p1.z() + u*(p2.z() - p1.z())
            );
            if (u > VSMALL && u < 1.0)
            { // closest point fall within the line segment
                minDist.insert(std::make_pair<scalar,point>(mag(p3 - px), px));
            }
        }

        if (minDist.empty())
        {
            Info<< "No nearest point found\n";
            // TODO use topo in order to find the next point
            return topo->getBoundaryPoint(pointRef);
        }
        else
        {
            return minDist.begin()->second;
        }
    }
    else
    {
        Info<< "Error feature edge point\n";
        return topo->getBoundaryPoint(pointRef);
    }
}

void Foam::featureEdgePoint::addPoint(const Foam::label &pointRef)
{    
    pointLinked_.insert(pointRef);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
