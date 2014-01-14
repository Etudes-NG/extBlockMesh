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

protected:

    //- Protected data

        //- Triangles faces connected to this point at beginning
        std::set<std::set<label> > triangles_;

        //- Triangles faces connected to this point during smoothing
        std::set<std::set<label> > trianglesNew_;

        // TODO remove storage of initial point
        //- Initial point
        point initialPoint_;


    //- Protected member functions

        //- Test triangle
        bool isOnTriangle
        (
            const label &refP2,
            const label &refP3,
            const point &p,
            point &out,
            const label &ref
        ) const;

        //- Change feature edge linked points and get optimal point
        virtual point changeFeatureEdgeLinkedsPoint
        (
            const label &newRef,
            const point &guessedPoint
        );

        //- Change boundary linked faces and get optimal point
        point changeBoundaryPointLinkedFaces
        (
            const label &newRef,
            const point &guessedPoint
        );

        //- Initial boundary point
        point getInitialPoint(const label &ref) const;

public:
    //- Constructors

        //- Construct from
        boundaryPoint
        (
            const std::set<std::set<Foam::label> > &triangles,
            const point &initialPoint,
            blockMeshTopology *topo
        );

    //- Destructor
        ~boundaryPoint();

    //- Member functions

        //- get optimal point with repect of point topo
        point smoothedPoint
        (
            const Foam::point &guessedPoint,
            const label &pointRef
        );

        //- Get map of neibor features edge points sorted by minimal distance
        std::map<scalar,point> mapNeiborFeaturePts
        (
            const point &guessedPoint,
            const label &pointRef
        );

        //- Get map of neibor boundary points sorted by minimal distance
        std::map<scalar,point> mapBoundaryFeaturePts
        (
            const point &guessedPoint,
            const label &pointRef
        );

        //- Get optimal boundary point from guessed
        point projectedBndPoint(const point &guessedPoint, const label &ref);


        //- Get optimal point from guessed
        virtual point getFeatureEdgePoint
        (
            const point &guessedPoint,
            const label &ref
        );

        virtual point getBoundaryPoint
        (
            const point &guessedPoint,
            const label &ref
        );

        //- Get triangles linked
        std::set<std::set<label> > getTrianglesLinked() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // BOUNDARYPOINT_H

// ************************************************************************* //

