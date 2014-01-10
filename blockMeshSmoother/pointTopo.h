#ifndef POINTTOPO_H
#define POINTTOPO_H

#include "blockMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockMeshTopology;

/*---------------------------------------------------------------------------*\
             Class pointTopo Declaration
\*---------------------------------------------------------------------------*/

class pointTopo
{
    //- Private data

        //-


    //- Private member functions


public:
    //- Constructors

        //- Construct from
        pointTopo();

    //- Destructor
        virtual ~pointTopo();

    //- Member functions

        //- get optimal point with repect of point topo
        virtual point smoothedPoint
        (
            const point &guessedPoint,
            const blockMesh *blocks,
            const label &pointRef,
            const blockMeshTopology *topo
        ) const {return guessedPoint;}

        virtual void addPoint(const label &addedPoint){}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // POINTTOPO_H

// ************************************************************************* //


