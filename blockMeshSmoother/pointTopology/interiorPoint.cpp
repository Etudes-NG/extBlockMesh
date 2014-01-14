#include "interiorPoint.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interiorPoint::interiorPoint(blockMeshTopology *topo)
:
    pointTopo(topo)
{
}

Foam::interiorPoint::~interiorPoint()
{

}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::point Foam::interiorPoint::smoothedPoint
(
    const Foam::point &guessedPoint,
    const Foam::label &pointRef
)
{
    // No constraints the new point is the guessed point
    return guessedPoint;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
