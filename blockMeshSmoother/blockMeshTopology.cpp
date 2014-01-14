#include "blockMeshTopology.h"

#include "featureEdgePoint.h"
#include "interiorPoint.h"
#include "boundaryPoint.h"
#include "cornerPoint.h"


#include <algorithm>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockMeshTopology::blockMeshTopology
(
    blockMesh *blocks,
    const dictionary &topoDict
)
    :
    pointTopo_(List<pointTopo*>(blockMeshPtr_->points().size())),
    featureAngle_
    (
      acos(-1.) - readScalar(topoDict.lookup("featureAngle"))/180*acos(-1.)
    ),
    minNbFeatureEdge_(readLabel(topoDict.lookup("minEdgeForFeature"))),
    pointsType_(List<pointType>(blockMeshPtr_->points().size(), internal)),
    blockMeshPtr_(blocks),
    boundaryFacePoints_(blocks->patches().size())
{
    boundaryFacePointsLabels();

    searchFeatureEdges(getBoundaryFacesPoints());
    initialiseBoundaryPoint();
}

Foam::blockMeshTopology::~blockMeshTopology()
{
    forAll (pointTopo_, ptI)
    {
        delete pointTopo_[ptI];
    }
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::blockMeshTopology::initialiseBoundaryPoint()
{
    const label nbPoints(blockMeshPtr_->points().size());
    List<std::set<std::set<label> > > pointTriangles(nbPoints);
    List<std::set<label> > linkedPointFace(nbPoints);
    std::set<label> boundaryPts;

    // Search triangles faces connected to point
    forAll (blockMeshPtr_->patches(), patchI)
    {
        forAll (blockMeshPtr_->patches()[patchI], faceI)
        {
            const labelList facePoints
            (
                blockMeshPtr_->patches()[patchI][faceI].pointsLabel()
            );
            forAll (facePoints, pointI)
            {
                // Search the connected triangle faces
                std::set<label> triangle1;
                triangle1.insert(facePoints[(pointI + 1) % 4]);
                triangle1.insert(facePoints[(pointI + 2) % 4]);

                std::set<label> triangle2;
                triangle2.insert(facePoints[(pointI + 2) % 4]);
                triangle2.insert(facePoints[(pointI + 3) % 4]);

                pointTriangles[facePoints[pointI]].insert(triangle1);
                pointTriangles[facePoints[pointI]].insert(triangle2);

                if (pointsType_[facePoints[pointI]] != featureEdge)
                { // Point is not a feature point

                    linkedPointFace[facePoints[pointI]].insert
                    (
                        facePoints[(pointI + 1) % 4]
                    );
                    linkedPointFace[facePoints[pointI]].insert
                    (
                        facePoints[(pointI + 3) % 4]
                    );
                    boundaryPts.insert(facePoints[pointI]);

                    pointsType_[facePoints[pointI]] = boundary;
                }
            }
        }
    }

    // Initialise list of pointer
    forAll (blockMeshPtr_->points(), ptI)
    {
        switch (pointsType_[ptI])
        {
        case featureEdge:
            pointTopo_[ptI] = new featureEdgePoint
            (
                featurePointConnections_[ptI],
                pointTriangles[ptI],
                blockMeshPtr_->points()[ptI],
                this
            );
            break;
        case internal:
            pointTopo_[ptI] = new interiorPoint
            (
                this
            );
            break;
        case corner:
            pointTopo_[ptI] = new cornerPoint
            (
                pointTriangles[ptI],
                blockMeshPtr_->points()[ptI],
                this
            );
            break;
        case boundary:
            pointTopo_[ptI] = new boundaryPoint
            (
                pointTriangles[ptI],
                blockMeshPtr_->points()[ptI],
                this
            );
            break;
        default:
            Info<< "Error wrong point topo!\n";
            break;
        }
    }
}


void Foam::blockMeshTopology::searchFeatureEdges
(
    const List<List<std::set<label> > > &bndFacesPoints
)
{
    const label nbPatch(blockMeshPtr_->patches().size());
    List<List<std::set<std::set<label> > > > bndFaceConnectEdges(nbPatch);
    List<List<std::set<std::set<label> > > > bndFaceNeiboor(nbPatch);

    // Search connected faces
    forAll (blockMeshPtr_->patches(), patchI)
    {
        const label nbFaces(blockMeshPtr_->patches()[patchI].size());

        bndFaceNeiboor[patchI].resize(nbFaces);
        bndFaceConnectEdges[patchI].resize(nbFaces);

        forAll (blockMeshPtr_->patches()[patchI], faceI)
        { // for all faces of this patch

            forAll (blockMeshPtr_->patches()[patchI], faceJ)
            { // Test all faces of this patch
                std::set<label> commonPts;
                std::set_intersection
                (
                    bndFacesPoints[patchI][faceI].begin(),
                    bndFacesPoints[patchI][faceI].end(),
                    bndFacesPoints[patchI][faceJ].begin(),
                    bndFacesPoints[patchI][faceJ].end(),
                    std::inserter(commonPts, commonPts.begin())
                );

                if (commonPts.size() == 2)
                { // Faces are connected
                    bndFaceNeiboor[patchI][faceI].insert
                    (
                        bndFacesPoints[patchI][faceJ]
                    );
                    bndFaceConnectEdges[patchI][faceI].insert(commonPts);

                    const scalar angle
                    (
                        connectedFaceAngles
                        (
                            patchI,
                            faceI,
                            faceJ,
                            *commonPts.begin(),
                            *commonPts.rbegin()
                        )
                    );

                    if (angle < featureAngle_)
                    { // New feature edge

                        insertFeaturePointConnection
                        (
                            *commonPts.begin(),
                            *commonPts.rbegin()
                        );
                        pointsType_[*commonPts.begin()] = featureEdge;
                        pointsType_[*commonPts.rbegin()] = featureEdge;
                    }
                }
            }
        }
    }

    // Store list of feature edges
    forAll (blockMeshPtr_->patches(), patchI)
    {
        forAll (blockMeshPtr_->patches()[patchI], faceI)
        {
            if(bndFaceNeiboor[patchI][faceI].size() < 4)
            { // One or more edge is a feature edge

                // Search edges not connected
                std::set<std::set<label> > faceEdge;
                const labelList ptLabel
                (
                    blockMeshPtr_->patches()[patchI][faceI].pointsLabel()
                );
                forAll (ptLabel, ptI)
                {
                    std::set<label> edge;
                    edge.insert(ptLabel[ptI]);
                    edge.insert(ptLabel[(ptI + 1) % 4]);
                    faceEdge.insert(edge);
                }

                // nb     edgeDef  pt
                std::set<std::set<label> > bndEdge;
                std::set_difference
                (
                    faceEdge.begin(),
                    faceEdge.end(),
                    bndFaceConnectEdges[patchI][faceI].begin(),
                    bndFaceConnectEdges[patchI][faceI].end(),
                    std::inserter(bndEdge, bndEdge.begin())
                );

                for // all feature edge of this face
                (
                    std::set<std::set<label> >::iterator
                        edgeI = bndEdge.begin();
                    edgeI != bndEdge.end();
                    ++edgeI
                )
                {
                    insertFeaturePointConnection
                    (
                        *(*edgeI).begin(),
                        *(*edgeI).rbegin()
                    );

                    pointsType_[*(*edgeI).begin()] = featureEdge;
                    pointsType_[*(*edgeI).rbegin()] = featureEdge;
                }
            }
        }
    }

    // FIXME add analyse of featureEdgeSet for nb of edge for each curve

    for
    (
        std::map<label, std::set<label> >::iterator ptI =
            featurePointConnections_.begin();
        ptI != featurePointConnections_.end();
        ++ptI
    )
    {
        if (ptI->second.size() != 2)
        { // Feature point is a corner
            pointsType_[ptI->first] = corner;
        }
    }
}

Foam::List<Foam::List<std::set<Foam::label> > >
Foam::blockMeshTopology::getBoundaryFacesPoints()
{
    const label nbPatch(blockMeshPtr_->patches().size());
    //patch  face  points labels
    List<List<std::set<label> > > bndFacesPoints(nbPatch);

    forAll (blockMeshPtr_->patches(), patchI)
    {
        const label nbFaces(blockMeshPtr_->patches()[patchI].size());
        bndFacesPoints[patchI].resize(nbFaces);

        forAll (blockMeshPtr_->patches()[patchI], faceI)
        {
            const labelList facePoints
            (
                blockMeshPtr_->patches()[patchI][faceI].pointsLabel()
            );
            std::set<label> facePts;
            forAll (facePoints, pointI)
            {
                facePts.insert(facePoints[pointI]);
            }

            bndFacesPoints[patchI][faceI] = facePts;
        }
    }
    return bndFacesPoints;
}

Foam::label Foam::blockMeshTopology::oppositePts
(
    const label &patch,
    const label &face,
    const label &contactPoints,
    const label &opposite
) const
{
    const label dist
    (
        std::find
        (
            boundaryFacePoints_[patch][face].begin(),
            boundaryFacePoints_[patch][face].end(),
            contactPoints
        ) -
        boundaryFacePoints_[patch][face].begin()
    );
    if(boundaryFacePoints_[patch][face][(dist + 1) % 4] == opposite)
    {
        return boundaryFacePoints_[patch][face][(dist + 3) % 4];
    }
    else
    {
        return boundaryFacePoints_[patch][face][(dist + 1) % 4];
    }
}

Foam::scalar Foam::blockMeshTopology::connectedFaceAngles
(
    const label &patch,
    const label &face1,
    const label &face2,
    const label &contactPoint1,
    const label &contactPoint2
) const
{
    const label A(oppositePts(patch, face1, contactPoint1, contactPoint2));
    const label B(oppositePts(patch, face2, contactPoint1, contactPoint2));

    const scalar alpha
    (
        (
            angleBetweenPlanes(A, B, contactPoint1, contactPoint2) +
            angleBetweenPlanes(A, B, contactPoint2, contactPoint1)
        )/2
    );

    return alpha;
}

Foam::scalar Foam::blockMeshTopology::angleBetweenPlanes
(
    const label &A,
    const label &B,
    const label &contactPoint1,
    const label &contactPoint2
) const
{
    const point v0B
    (
        blockMeshPtr_->points()[B] - blockMeshPtr_->points()[contactPoint1]
    );
    const point v01
    (
        blockMeshPtr_->points()[contactPoint2] -
        blockMeshPtr_->points()[contactPoint1]
    );
    const point v0A
    (
        blockMeshPtr_->points()[A] - blockMeshPtr_->points()[contactPoint1]
    );

    const point nt1(v0B ^ v01);
    const point nt2(v0A ^ v01);

    const scalar Acos((nt1 & nt2)/(mag(nt1)*mag(nt2) + VSMALL));
    if (Acos > 1.0)
    {
        return acos(1.0);
    }
    else if (Acos < -1.0)
    {
        return acos(-1.0);
    }
    else
    {
        return acos(Acos);
    }
}

void Foam::blockMeshTopology::insertFeaturePointConnection
(
    const Foam::label &pt1,
    const Foam::label &pt2
)
{
    std::map<label, std::set<label> >::iterator match =
            featurePointConnections_.find(pt1);

    if (match != featurePointConnections_.end())
    {
        match->second.insert(pt2);
    }
    else
    {
        std::set<label> lk;
        lk.insert(pt2);
        featurePointConnections_[pt1] = lk;
    }

    match = featurePointConnections_.find(pt2);

    if (match != featurePointConnections_.end())
    {
        match->second.insert(pt1);
    }
    else
    {
        std::set<label> lk;
        lk.insert(pt1);
        featurePointConnections_[pt2] = lk;
    }
}

void Foam::blockMeshTopology::boundaryFacePointsLabels()
{
    forAll (blockMeshPtr_->patches(), patchI)
    {
        boundaryFacePoints_[patchI].resize
        (
            blockMeshPtr_->patches()[patchI].size()
        );
        forAll (blockMeshPtr_->patches()[patchI], faceI)
        {

            boundaryFacePoints_[patchI][faceI] =
                blockMeshPtr_->patches()[patchI][faceI].pointsLabel();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::point Foam::blockMeshTopology::optimalPoint
(
    const label &pointRef,
    const point &guessedPoint
)
{
    return pointTopo_[pointRef]->smoothedPoint(guessedPoint, pointRef);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //

