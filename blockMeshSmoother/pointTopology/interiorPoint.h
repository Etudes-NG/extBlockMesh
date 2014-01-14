#ifndef INTERIORPOINT_H
#define INTERIORPOINT_H

#include "pointTopo.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class pointTopo Declaration
\*---------------------------------------------------------------------------*/

class interiorPoint : public pointTopo
{
    //- Private data

        //-


    //- Private member functions


public:
    //- Constructors

        //- Construct from
        interiorPoint(blockMeshTopology *topo);

    //- Destructor
        ~interiorPoint();

    //- Member functions

        //- Point smoothed with respect of topology constraints
        point smoothedPoint
        (
            const point &guessedPoint,
            const label &pointRef
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // INTERIORPOINT_H

// ************************************************************************* //

