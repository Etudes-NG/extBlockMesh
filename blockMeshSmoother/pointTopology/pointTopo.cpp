#include "pointTopo.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointTopo::pointTopo(blockMeshTopology *topo)
    :
      topo_(topo)
{
}

Foam::pointTopo::~pointTopo()
{

}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::set<Foam::label> Foam::pointTopo::getPointLinked() const
{
    FatalErrorIn("getPointLinked()")
        << "Accessed from a non feature edge point\n"
        << nl
        << exit(FatalError);

    return std::set<label>();
}

Foam::point Foam::pointTopo::smoothedPoint
(
    const Foam::point &guessedPoint,
    const Foam::label &pointRef
)
{
    FatalErrorIn("smoothedPoint()")
        << "Accessed from a non feature edge point\n"
        << nl
        << exit(FatalError);

    return point();
}

std::map<Foam::scalar, Foam::point> Foam::pointTopo::mapNeiborFeaturePts
(
    const Foam::point &guessedPoint,
    const Foam::label &pointRef
) const
{
    FatalErrorIn("mapNeiborFeaturePts()")
        << "Accessed from a non feature edge point\n"
        << nl
        << exit(FatalError);


    return std::map<scalar,point>();
}

std::set<std::set<Foam::label> > Foam::pointTopo::getTrianglesLinked() const
{
    FatalErrorIn("getTrianglesLinked()")
        << "Accessed from a non feature edge point\n"
        << nl
        << exit(FatalError);

    return std::set<std::set<Foam::label> >();
}

Foam::point Foam::pointTopo::getFeatureEdgePoint
(
    const Foam::point &guessedPoint,
    const Foam::label &ref
)
{
    FatalErrorIn("getFeatureEdgePoint()")
        << "Accessed from a non feature edge point\n"
        << nl
        << exit(FatalError);

    return point();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
