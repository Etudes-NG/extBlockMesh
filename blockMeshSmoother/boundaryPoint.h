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

        //- Point topology
        std::set<std::set<label> > triangles_;

        //- Normal point
        std::set<label> normalPoint_;

    //- Private member functions

public:
    //- Constructors

        //- Construct from
        boundaryPoint
        (
            const std::set<std::set<Foam::label> > &triangles,
            const std::set<label> &normalPoint
        );

    //- Destructor
        ~boundaryPoint();

    //- Member functions

        //- get optimal point with repect of point topo
        point smoothedPoint
        (
            const Foam::point &guessedPoint,
            const blockMesh *blocks,
            const label &pointRef,
            const blockMeshTopology *topo
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // BOUNDARYPOINT_H

// ************************************************************************* //

