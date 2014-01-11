#include "boundaryPoint.h"

#include "blockMeshTopology.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryPoint::boundaryPoint
(const std::set<std::set<Foam::label> > &triangles)
    :
      triangles_(triangles)
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
    const label &pointRef
)
{
    // TODO: compute guessed point from normal and topo
    return guessedPoint;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
