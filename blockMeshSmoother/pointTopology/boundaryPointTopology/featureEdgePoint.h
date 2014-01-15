#ifndef FEATUREEDGEPOINT_H
#define FEATUREEDGEPOINT_H

#include "boundaryPoint.h"

#include <set>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class pointTopo Declaration
\*---------------------------------------------------------------------------*/

class featureEdgePoint : public boundaryPoint
{
    //- Private data

        //- Point link at beginning
        std::set<label> pointLinked_;

        //- Point link during move
        std::set<label> pntLinkedNew_;

        //- Initial point label
        label featurePtLabel_;

        //- Point label during smoothing
        label newFeaturePtLabel_;

        // TODO change "new" by pointer of pointTopo

    //- Private member functions

        //- Get linked point of feature edge point
        std::set<label> getPointLinked() const;

public:
    //- Constructors

        //- Construct from
        featureEdgePoint
        (
            const std::set<label> &pointLinked,
            const std::set<std::set<Foam::label> > &triangles,
            const point &initialPoint,
            const label &initialLabel,
            blockMeshTopology *topo
        );

    //- Destructor
        ~featureEdgePoint();

    //- Member functions

        //- Get optimal point with repect of point topo
        point smoothedPoint(const point &guessedPoint);

        //- Get map of points sorted by minimal distance
        std::map<scalar,point> mapNeiborFeaturePts
        (
            const point &guessedPoint
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

#endif // FEATUREEDGEPOINT_H

// ************************************************************************* //

