#ifndef POINTTOPO_H
#define POINTTOPO_H

#include "blockMesh.H"

#include <map>
#include <set>

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
        virtual ~pointTopo() = 0;

    //- Member functions

        //- get optimal point with repect of point topo
        virtual point smoothedPoint
        (
            const point &guessedPoint,
            const label &pointRef
        )
        {
            Info<< "call wrong point topo\n";
            return point();
        }

        virtual std::map<scalar,point> minDist
        (
            const point &guessedPoint,
            const label &pointRef
        ) const {return std::map<scalar,point>();}

        //- Get linked point of feature edge point
        virtual std::set<label> getPointLinked() const
        {return std::set<label>();}

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // POINTTOPO_H

// ************************************************************************* //


