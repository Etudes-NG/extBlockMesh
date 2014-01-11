#include "blockMeshTopology.h"

#include "featureEdgePoint.h"
#include "interiorPoint.h"
#include "boundaryPoint.h"

#include <algorithm>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockMeshTopology::blockMeshTopology
(
    const blockMesh *blocks,
    const dictionary &topoDict
)
    :
      pointTopo_(List<pointTopo*>(blocks->points().size())),
      featureAngle_
      (
          acos(-1.) - readScalar(topoDict.lookup("featureAngle"))/180*acos(-1.)
      ),
      minNbFeatureEdge_(readLabel(topoDict.lookup("minEdgeForFeature")))
{
    searchFeatureEdges(getBoundaryFacesPoints(blocks), blocks);
    initialiseBoundaryPoint(blocks);

    Info<< "There is " << label(featureEdgeSet_.size()) << "feature edges" << nl;
    Info<< "There is " << label(featurePts_.size()) << "feature points" << nl;
}

Foam::blockMeshTopology::~blockMeshTopology()
{
    forAll (pointTopo_, ptI)
    {
        delete pointTopo_[ptI];
    }
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::blockMeshTopology::initialiseBoundaryPoint(const blockMesh *blocks)
{
    List<std::set<std::set<label> > > pointTriangles(blocks->points().size());
    List<std::set<label> > linkedPointFace(blocks->points().size());
    std::set<label> boundaryPts;

    // Search triangles connected to point
    forAll (blocks->patches(), patchI)
    {
        forAll (blocks->patches()[patchI], faceI)
        {
            const labelList facePoints
            (
                blocks->patches()[patchI][faceI].pointsLabel()
            );
            forAll (facePoints, pointI)
            {
                if (featurePts_.find(facePoints[pointI]) == featurePts_.end())
                { // Point is not a feature point

                    std::set<label> triangle1;
                    triangle1.insert(facePoints[(pointI + 1) % 4]);
                    triangle1.insert(facePoints[(pointI + 2) % 4]);

                    std::set<label> triangle2;
                    triangle2.insert(facePoints[(pointI + 2) % 4]);
                    triangle2.insert(facePoints[(pointI + 3) % 4]);

                    pointTriangles[facePoints[pointI]].insert(triangle1);
                    pointTriangles[facePoints[pointI]].insert(triangle2);

                    linkedPointFace[facePoints[pointI]].insert
                    (
                        facePoints[(pointI + 1) % 4]
                    );
                    linkedPointFace[facePoints[pointI]].insert
                    (
                        facePoints[(pointI + 3) % 4]
                    );
                    boundaryPts.insert(facePoints[pointI]);
                }
            }
        }
    }

    // Initialise list of pointer
    forAll (blocks->points(), ptI)
    {
        if (boundaryPts.find(ptI) != boundaryPts.end())
        { // Point is a boundaryPoint, initialise it
            pointTopo_[ptI] = new boundaryPoint(pointTriangles[ptI]);
        }
        else if (featurePts_.find(ptI) != featurePts_.end())
        { // Point is a feature point

            pointTopo_[ptI] = new featureEdgePoint
            (
                featurePointConnections_[ptI],
                this
            );
        }
        else
        { // Point is not a boundary point and not a feature point
            pointTopo_[ptI] = new interiorPoint();
        }
    }
}


void Foam::blockMeshTopology::searchFeatureEdges
(
    const List<List<std::set<label> > > &bndFacesPoints,
    const blockMesh *blocks
)
{
    const label nbPatch(blocks->patches().size());
    List<List<std::set<std::set<label> > > > bndFaceConnectEdges(nbPatch);
    List<List<std::set<std::set<label> > > > bndFaceNeiboor(nbPatch);

    // Search connected faces
    forAll (blocks->patches(), patchI)
    {
        const label nbFaces(blocks->patches()[patchI].size());

        bndFaceNeiboor[patchI].resize(nbFaces);
        bndFaceConnectEdges[patchI].resize(nbFaces);

        forAll (blocks->patches()[patchI], faceI)
        { // for all faces of this patch

            forAll (blocks->patches()[patchI], faceJ)
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
                            *commonPts.rbegin(),
                            blocks
                        )
                    );

                    if (angle < featureAngle_)
                    { // New feature edge

                        insertFeaturePointConnection
                        (
                            *commonPts.begin(),
                            *commonPts.rbegin()
                        );
                        featurePts_.insert(*commonPts.begin());
                        featurePts_.insert(*commonPts.rbegin());

                        featureEdgeSet_.insert(commonPts);
                    }
                }
            }
        }
    }

    // Store list of feature edges
    forAll (blocks->patches(), patchI)
    {
        forAll (blocks->patches()[patchI], faceI)
        {
            if(bndFaceNeiboor[patchI][faceI].size() < 4)
            { // One or more edge is a feature edge

                // Search edges not connected
                std::set<std::set<label> > faceEdge;
                const labelList ptLabel
                (
                    blocks->patches()[patchI][faceI].pointsLabel()
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

                for
                (
                    std::set<std::set<label> >::iterator
                        iter = bndEdge.begin();
                    iter != bndEdge.end();
                    ++iter
                )
                {
                    insertFeaturePointConnection
                    (
                        *(*iter).begin(),
                        *(*iter).rbegin()
                    );
                    featurePts_.insert(*(*iter).begin());
                    featurePts_.insert(*(*iter).rbegin());
                    featureEdgeSet_.insert(*iter);
                }
            }
        }
    }
}

Foam::List<Foam::List<std::set<Foam::label> > >
Foam::blockMeshTopology::getBoundaryFacesPoints(const blockMesh *blocks)
{
    const label nbPatch(blocks->patches().size());
    //patch  face  points labels
    List<List<std::set<label> > > bndFacesPoints(nbPatch);

    forAll (blocks->patches(), patchI)
    {
        const label nbFaces(blocks->patches()[patchI].size());
        bndFacesPoints[patchI].resize(nbFaces);

        forAll (blocks->patches()[patchI], faceI)
        {
            const labelList facePoints
            (
                blocks->patches()[patchI][faceI].pointsLabel()
            );
            std::set<label> facePts;
            forAll (facePoints, pointI)
            {
                facePts.insert(facePoints[pointI]);

                // Store boundary point
                bndPoints_[facePoints[pointI]] =
                        blocks->points()[facePoints[pointI]];
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
    const label &opposite,
    const blockMesh *blocks
) const
{
    const labelList ptLabel(blocks->patches()[patch][face].pointsLabel());

    const label dist
    (
        std::find(ptLabel.begin(), ptLabel.end(), contactPoints) -
        ptLabel.begin()
    );
    if(ptLabel[(dist + 1) % 4] == opposite)
    {
        return ptLabel[(dist + 3) % 4];
    }
    else
    {
        return ptLabel[(dist + 1) % 4];
    }
}

Foam::scalar Foam::blockMeshTopology::connectedFaceAngles
(
    const label &patch,
    const label &face1,
    const label &face2,
    const label &contactPoint1,
    const label &contactPoint2,
    const blockMesh *blocks
) const
{
    const label A
    (
        oppositePts(patch, face1, contactPoint1, contactPoint2, blocks)
    );
    const label B
    (
        oppositePts(patch, face2, contactPoint1, contactPoint2, blocks)
    );

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
    const Foam::label &A,
    const Foam::label &B,
    const Foam::label &contactPoint1,
    const Foam::label &contactPoint2
) const
{
    const point v0B(bndPoints_.at(B) - bndPoints_.at(contactPoint1));
    const point v01(bndPoints_.at(contactPoint2) -bndPoints_.at(contactPoint1));
    const point v0A(bndPoints_.at(A) - bndPoints_.at(contactPoint1));

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

