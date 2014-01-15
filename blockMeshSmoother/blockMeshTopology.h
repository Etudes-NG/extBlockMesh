#ifndef BLOCKMESHTOPOLOGY_H
#define BLOCKMESHTOPOLOGY_H

#include "blockMesh.H"

#include "pointTopo.h"

#include <set>
#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class blockMeshTopology Declaration
\*---------------------------------------------------------------------------*/

class blockMeshTopology
{
    // Private data

        //- Type of points
        enum pointType
        {
            internal, featureEdge, corner, boundary
        };

        //- Point topology
        List<pointTopo*> pointTopo_;

        //- angle target for feature edge
        scalar featureAngle_;

        //- Min feature edge lenght
        label minNbFeatureEdge_;

        //- Features point connections
        std::map<label, std::set<label> > featurePointConnections_;

        //- List of points with type
        List<pointType> pointsType_;

        //- Pointer of blockMesh
        blockMesh *blockMeshPtr_;

        //- TODO add storage of cell point conectivity (from blockMeshSmoother)

        //- Face points label in list
        labelListListList boundaryFacePoints_;


    //- Private member functions

        //- Add boudary point triangles
        void initialiseBoundaryPoint();

        //- Get patch faces neiboor
        void searchFeatureEdges
        (
            const List<List<std::set<label> > > &bndFacesPoints
        );

        //- Search boundary faces edges
        List<List<std::set<label> > > getBoundaryFacesPoints();

        //- Search point normal of boundary faces
        List<std::set<label> > pointLinks();

        //- Get normal point of boundary faces
        label oppositePts
        (
            const label &patch,
            const label &face,
            const label &contactPoints,
            const label &opposite
        ) const;

        scalar connectedFaceAngles
        (
            const label &patch,
            const label &face1,
            const label &face2,
            const label &contactPoint1,
            const label &contactPoint2
         ) const;

        scalar angleBetweenPlanes
        (
            const label &A,
            const label &B,
            const label &contactPoint1,
            const label &contactPoint2
        ) const;

        //- Insert feature point connection
        void insertFeaturePointConnection(const label &pt1, const label &pt2);

        //- Store Face points label in list
        void boundaryFacePointsLabels();
public:
    //- Constructors

        //- Construct from
        blockMeshTopology(blockMesh *blocks, const dictionary &topoDict);

    //- Destructor
        ~blockMeshTopology();


    //- Member functions

        //- Get optimal point from gessed and normal point
        point optimalPoint
        (
            const label &pointRef,
            const point &guessedPoint
        );

        //- Get boundary point data
        point &getBndPointCoord(const label &ref) const
        {
            return pointTopo_[ref]->getboundaryPoint();
        }

        pointTopo *getPointTopoPtr(const label &ref) const
        {
            return pointTopo_[ref];
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // BLOCKMESHTOPOLOGY_H

// ************************************************************************* //
