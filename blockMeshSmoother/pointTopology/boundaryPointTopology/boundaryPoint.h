#ifndef BOUNDARYPOINT_H
#define BOUNDARYPOINT_H

#include "pointTopo.h"

#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class boundaryPoint Declaration
\*---------------------------------------------------------------------------*/

class boundaryPoint : public pointTopo
{
    //- Private data

    //- Private member functions

        //- Test triangle
        bool isOnTriangle
        (
            const label &refP2,
            const label &refP3,
            const point &p,
            point &out,
            const label &ref
       ) const;

        //- Get triangles linked
        std::set<std::set<label> > getTrianglesLinked() const;

        //- Change linked points and get optimal point
        point changeLinkedsPoint
        (
            const label &newRef,
            const point &guessedPoint
        );

public:
    //- Constructors

        //- Construct from
        boundaryPoint
        (
            const std::set<std::set<Foam::label> > &triangles,
            blockMeshTopology *topo
        );

    //- Destructor
        ~boundaryPoint();

    //- Member functions

        //- get optimal point with repect of point topo
        point smoothedPoint
        (
            const Foam::point &guessedPoint,
            const label &pointRef
        );

        //- Get map of points sorted by minimal distance
        std::map<scalar,point> minDist
        (
            const point &guessedPoint,
            const label &pointRef
        );

        //- Get optimal boundary point from guessed
        point projectedBndPoint(const point &guessedPoint, const label &ref);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // BOUNDARYPOINT_H

// ************************************************************************* //

