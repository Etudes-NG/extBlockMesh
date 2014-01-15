#ifndef CORNERPOINT_H
#define CORNERPOINT_H

#include "boundaryPoint.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cornerPoint Declaration
\*---------------------------------------------------------------------------*/


class cornerPoint : public boundaryPoint
{
public:
    //- Constructors

        //- Construct from
        cornerPoint
        (
            const std::set<std::set<Foam::label> > &triangles,
            const point &initialPoint,
            const label &initialLabel,
            blockMeshTopology *topo
        );

    //- Destructor
        ~cornerPoint();

    //- Member functions

        //- get optimal point with repect of point topo
        virtual point smoothedPoint
        (
            const point &guessedPoint,
            const label &pointRef
        ) const;

        //- Change feature edge linked points and get optimal point
        point changeFeatureEdgeLinkedsPoint
        (
            const label &newRef,
            const point &guessedPoint
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // CORNERPOINT_H

// ************************************************************************* //

