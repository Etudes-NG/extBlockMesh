#include "cornerPoint.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cornerPoint::cornerPoint
(
    const std::set<std::set<Foam::label> > &triangles,
    const point &initialPoint,
    const label &initialLabel,
    blockMeshTopology *topo
)
    :
      boundaryPoint(triangles, initialPoint, initialLabel, topo)
{
}

Foam::cornerPoint::~cornerPoint()
{

}


// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::cornerPoint::smoothedPoint
(
    const point &guessedPoint,
    const label &pointRef
) const
{
    // No change as the point is fixed
    return initialPoint_;
}

Foam::point Foam::cornerPoint::changeFeatureEdgeLinkedsPoint
(
    const Foam::label &newRef,
    const Foam::point &guessedPoint
)
{
    FatalErrorIn("changeFeatureEdgeLinkedsPoint(guessedPoint, ref)")
        << "Accessed from a non feature edge point\n"
        << nl
        << exit(FatalError);

    return point();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //

