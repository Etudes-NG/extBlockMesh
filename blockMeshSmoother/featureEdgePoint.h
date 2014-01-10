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

        //- Point link 1
        std::set<label> pointLinked_;


    //- Private member functions


public:
    //- Constructors

        //- Construct from
        featureEdgePoint();

    //- Destructor
        ~featureEdgePoint();

    //- Member functions

        //- Start smoothing
        point smoothedPoint
        (
            const Foam::point &guessedPoint,
            const blockMesh *blocks,
            const label &pointRef,
            const blockMeshTopology *topo
        ) const;

        void addPoint(const label &pointRef);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // FEATUREEDGEPOINT_H

// ************************************************************************* //

