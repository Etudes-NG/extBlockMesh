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
        interiorPoint();

    //- Destructor
        ~interiorPoint();

    //- Member functions

        //- Start smoothing
        point smoothedPoint
        (
            const point &guessedPoint,
            const label &pointRef
        )
        {
            return guessedPoint;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // INTERIORPOINT_H

// ************************************************************************* //

