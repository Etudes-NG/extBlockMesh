#ifndef FEATUREEDGEPOINT_H
#define FEATUREEDGEPOINT_H

#include "pointTopo.h"

#include <set>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class pointTopo Declaration
\*---------------------------------------------------------------------------*/

class featureEdgePoint : public pointTopo
{
    //- Private data

        //- Point link at beginning
        std::set<label> pointLinked_;

        //- Point link when moving
        std::set<label> pointLinkedNew_;

        //- Pointer of topology
        blockMeshTopology *topo_;

    //- Private member functions

        //- Get linked point of feature edge point
        std::set<label> getPointLinked() const;

        //- Get optimal point from guessed
        point getPoint(const point &guessedPoint, const label &ref);

        //- Change linked points and get optimal point
        point changeLinkedsPoint
        (
            const label &newRef,
            const point &guessedPoint
        );

public:
    //- Constructors

        //- Construct from
        featureEdgePoint
        (
            const std::set<label> &pointLinked,
            blockMeshTopology *topo
        );

    //- Destructor
        ~featureEdgePoint();

    //- Member functions

        //- Start smoothing
        point smoothedPoint
        (
            const point &guessedPoint,
            const label &pointRef
        );

        std::map<scalar,point> minDist
        (
            const point &guessedPoint,
            const label &pointRef
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // FEATUREEDGEPOINT_H

// ************************************************************************* //

