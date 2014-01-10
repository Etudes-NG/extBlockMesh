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
    searchFeatureEdges(boundaryFaceEdge(blocks), blocks);
    initialiseBoundaryPoint(pointLinks(blocks), blocks);
}

Foam::blockMeshTopology::~blockMeshTopology()
{
    forAll (pointTopo_, ptI)
    {
        delete pointTopo_[ptI];
    }
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::blockMeshTopology::initialiseBoundaryPoint
(
    const List<std::set<label> > &pointsLink,
    const blockMesh *blocks
)
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
    Info<< "Nb boundary points = " << label(boundaryPts.size()) << endl;
    Info<< "Nb features points = " << label(featurePts_.size()) << endl;

    forAll (blocks->points(), ptI)
    {
        if (boundaryPts.find(ptI) != boundaryPts.end())
        { // Point is a boundaryPoint, initialise it


            // search the normal point
            std::set<label> normalPoint;
            std::set_difference
            (
                pointsLink[ptI].begin(),
                pointsLink[ptI].end(),
                linkedPointFace[ptI].begin(),
                linkedPointFace[ptI].end(),
                std::inserter(normalPoint, normalPoint.begin())
            );

            pointTopo_[ptI] = new Foam::boundaryPoint
            (
                pointTriangles[ptI],
                normalPoint
            );
        }
        else if (featurePts_.find(ptI) == featurePts_.end())
        { // Point is not a boundary point and not a feature point
            pointTopo_[ptI] = new interiorPoint();
        }
    }
}


Foam::labelListList Foam::blockMeshTopology::orderEdgesCurve
(
    std::set<std::set<label> > &featEdgSet
)
{
    labelListList featEdge;
    std::set<label> edgePoint;

    while (!featEdgSet.empty())
    { // While there is featureEdge not connected

        // create a new feature edge
        featEdge.append(labelList());

        // Add point 1 of edge 1 in feature edge def
        //                          patch    edge     pt0
        const label pt0(*featEdgSet.begin()->begin());
        featEdge[featEdge.size() - 1].append(pt0);
        const label pt1(*featEdgSet.begin()->rbegin());
        featEdge[featEdge.size() - 1].append(pt1);
        edgePoint.insert(pt0);
        edgePoint.insert(pt1);

        // erase edgeDef 1 from patch edges
        featEdgSet.erase(featEdgSet.begin());

        for
        (
            std::set<std::set<label> >::iterator edgeI = featEdgSet.begin();
            edgeI != featEdgSet.end();
            ++edgeI
        )
        { // for each edge

            for
            (
                std::set<std::set<label> >::iterator edgeK = featEdgSet.begin();
                edgeK != featEdgSet.end();
                ++edgeK
            )
            { // Test all edges
                const label fn(featEdge.size() - 1);
                if (*edgeK->begin() == *featEdge[fn].rbegin())
                {
                    const label ptn(*edgeK->rbegin());
                    featEdge[featEdge.size() - 1].append(ptn);
                    edgePoint.insert(ptn);
                    featEdgSet.erase(edgeK);
                }
                else if (*edgeK->rbegin() == *featEdge[fn].rbegin())
                {
                    const label ptn(*edgeK->begin());
                    featEdge[featEdge.size() - 1].append(ptn);
                    edgePoint.insert(ptn);
                    featEdgSet.erase(edgeK);
                }
            }
        }
    }

    return featEdge;
}


void Foam::blockMeshTopology::searchFeatureEdges
(
    const List<List<std::set<label> > > &bndFacesName,
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
                    bndFacesName[patchI][faceI].begin(),
                    bndFacesName[patchI][faceI].end(),
                    bndFacesName[patchI][faceJ].begin(),
                    bndFacesName[patchI][faceJ].end(),
                    std::inserter(commonPts, commonPts.begin())
                );

                if (commonPts.size() == 2)
                { // Faces share an edge
                    bndFaceNeiboor[patchI][faceI].insert
                    (
                        bndFacesName[patchI][faceJ]
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

    // Initialise list of pointer for feature point
    for
    (
        std::set<label>::iterator iter = featurePts_.begin();
        iter != featurePts_.end();
        ++iter
    )
    {
        pointTopo_[*iter] = new featureEdgePoint();
    }
}

Foam::List<Foam::List<std::set<Foam::label> > >
Foam::blockMeshTopology::boundaryFaceEdge(const blockMesh *blocks)
{
    const label nbPatch(blocks->patches().size());
    //patch  face  points labels
    List<List<std::set<label> > > bndFacesName(nbPatch);

    forAll (blocks->patches(), patchI)
    {
        const label nbFaces(blocks->patches()[patchI].size());
        bndFacesName[patchI].resize(nbFaces);

        std::set<std::set<label> > facesPointsSet;
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
                bndPoints_.insert
                (
                    std::make_pair<label, point>
                    (
                        facePoints[pointI],
                        blocks->points()[facePoints[pointI]]
                    )
                );
            }
            facesPointsSet.insert(facePts);

            bndFacesName[patchI][faceI] = facePts;
        }
    }
    return bndFacesName;
}

Foam::List<std::set<Foam::label> > Foam::blockMeshTopology::pointLinks
(
    const blockMesh *blocks
)
{
    List<std::set<label> > pointsLinks(blocks->points().size());
    forAll (blocks->cells(), cellI)
    {
        const labelList ptsLabels(blocks->cells()[cellI].pointsLabel());

        pointsLinks[ptsLabels[0]].insert(ptsLabels[1]);
        pointsLinks[ptsLabels[0]].insert(ptsLabels[3]);
        pointsLinks[ptsLabels[0]].insert(ptsLabels[4]);

        pointsLinks[ptsLabels[1]].insert(ptsLabels[0]);
        pointsLinks[ptsLabels[1]].insert(ptsLabels[2]);
        pointsLinks[ptsLabels[1]].insert(ptsLabels[5]);

        pointsLinks[ptsLabels[2]].insert(ptsLabels[1]);
        pointsLinks[ptsLabels[2]].insert(ptsLabels[6]);
        pointsLinks[ptsLabels[2]].insert(ptsLabels[3]);

        pointsLinks[ptsLabels[3]].insert(ptsLabels[2]);
        pointsLinks[ptsLabels[3]].insert(ptsLabels[7]);
        pointsLinks[ptsLabels[3]].insert(ptsLabels[0]);

        pointsLinks[ptsLabels[4]].insert(ptsLabels[7]);
        pointsLinks[ptsLabels[4]].insert(ptsLabels[5]);
        pointsLinks[ptsLabels[4]].insert(ptsLabels[0]);

        pointsLinks[ptsLabels[5]].insert(ptsLabels[4]);
        pointsLinks[ptsLabels[5]].insert(ptsLabels[6]);
        pointsLinks[ptsLabels[5]].insert(ptsLabels[1]);

        pointsLinks[ptsLabels[6]].insert(ptsLabels[5]);
        pointsLinks[ptsLabels[6]].insert(ptsLabels[7]);
        pointsLinks[ptsLabels[6]].insert(ptsLabels[2]);

        pointsLinks[ptsLabels[7]].insert(ptsLabels[6]);
        pointsLinks[ptsLabels[7]].insert(ptsLabels[4]);
        pointsLinks[ptsLabels[7]].insert(ptsLabels[3]);
    }

    return pointsLinks;
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
    const point &guessedPoint,
    const blockMesh *blocks
)
{
    return pointTopo_[pointRef]->smoothedPoint
    (
        guessedPoint,
        blocks,
        pointRef,
        this
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //

