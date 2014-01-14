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

    //- Private member functions

protected:
        //- Pointer of boundary points topology
        blockMeshTopology *topo_;

    //- Protected member functions


public:
    //- Constructors

        //- Construct from
        pointTopo(blockMeshTopology *topo);

    //- Destructor
        virtual ~pointTopo() = 0;

    //- Member functions

        //- Point smoothed with respect of topology constraints
        virtual point smoothedPoint
        (
            const point &guessedPoint,
            const label &pointRef
        );

        virtual std::map<scalar,point> mapNeiborFeaturePts
        (
            const point &guessedPoint,
            const label &pointRef
        ) const;

        //- Get linked point of feature edge point
        virtual std::set<label> getPointLinked() const;

        //- Get linked triangles faces of feature edge point
        virtual std::set<std::set<label> > getTrianglesLinked() const;


        //- Get optimal point from guessed
        virtual point getFeatureEdgePoint
        (
            const point &guessedPoint,
            const label &ref
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // POINTTOPO_H

// ************************************************************************* //


