#include "boundaryPoint.h"

#include "blockMeshTopology.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryPoint::boundaryPoint
(
    const std::set<std::set<Foam::label> > &triangles,
    const std::set<label> &normalPoint
)
    :
      triangles_(triangles),
      normalPoint_(normalPoint)
{
}

Foam::boundaryPoint::~boundaryPoint()
{
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::boundaryPoint::smoothedPoint
(
    const Foam::point &guessedPoint,
    const blockMesh *blocks,
    const label &pointRef,
    const blockMeshTopology *topo
) const
{
    // TODO: compute guessed point from normal and topo
    return guessedPoint;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
