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

        //- Point topology
        List<pointTopo*> pointTopo_;

        //- initial point data
        std::map<label,point> bndPoints_;

        //- angle target for feature edge
        scalar featureAngle_;

        //- Min feature edge lenght
        label minNbFeatureEdge_;

        //- Features point connections
        std::map<label, std::set<label> > featurePointConnections_;

        //- TODO add storage of cell point conectivity (from blockMeshSmoother)
        //- TODO store face points label in list


    //- Private member functions

        //- Initialisee the list for feature edge point
        std::set<label> initiliseFeatureEdgePoint
        (
            const labelListList &featureEdgeDef
        );

        //- Add boudary point triangles
        void initialiseBoundaryPoint
        (
            const std::set<label> &featurePts,
            const List<std::set<label> > &pointsLink,
            const Foam::blockMesh *blocks
        );

        //- Order edges forming a curve
        labelListList orderEdgesCurve(std::set<std::set<label> > &featEdgSet);

        //- Extract boundary of patchs
        void extractPatchBoundary
        (
            const List<List<std::set<std::set<label> > > > &bndFaceNeiboor,
            const List<List<std::set<std::set<label> > > > &bndFaceConnectEdges,
            const blockMesh *blocks
        );

        //- Get patch faces neiboor
        List<List<std::set<std::set<label> > > > getPatchFacesNeiboor
        (
            List<List<std::set<std::set<label> > > > &bndFaceNeiboor,
            const List<List<std::set<label> > > &bndFacesName,
            const blockMesh *blocks
        );

        //- Search boundary faces edges
        List<List<std::set<label> > > boundaryFacesEdges
        (
            List<std::set<std::set<label> > > &patchFacesPoints,
            const blockMesh *blocks
        );

        //- Search point normal of boundary faces
        List<std::set<label> > pointLinks(const blockMesh *blocks);

        //- Get normal point of boundary faces
        label oppositePts
        (
            const label &patch,
            const label &face,
            const label &contactPoints,
            const label &opposite,
            const blockMesh *blocks
        ) const;

        scalar connectedFaceAngles
        (
            const label &patch,
            const label &face1,
            const label &face2,
            const label &contactPoint1,
            const label &contactPoint2,
            const blockMesh *blocks
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
public:
    //- Constructors

        //- Construct from
        blockMeshTopology(const blockMesh *blocks, const dictionary &topoDict);

    //- Destructor
        ~blockMeshTopology();


    //- Member functions

        //- Get optimal point from gessed and normal point
        point optimalPoint
        (
            const label &pointRef,
            const point &guessedPoint,
            const blockMesh *blocks
        );

        //- Get boundary point data
        point getBoundaryPoint(const label &ref) const
        {
            return bndPoints_.at(ref);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // BLOCKMESHTOPOLOGY_H

// ************************************************************************* //
