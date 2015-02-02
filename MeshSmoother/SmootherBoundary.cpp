/*---------------------------------------------------------------------------*\
  extBlockMesh
  Copyright (C) 2014 Etudes-NG
  ---------------------------------
License
    This file is part of extBlockMesh.

    extBlockMesh is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    extBlockMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with extBlockMesh.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "SmootherBoundary.h"

#include "SmootherVertex.h"
#include "SmootherEdge.h"
#include "SmootherSurface.h"

#include "dictionary.H"
#include "polyMesh.H"
#include "Time.H"
#include <OFstream.H>
#include "unitConversion.H"

#include <sstream>

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::SmootherBoundary::analyseDict(dictionary &snapDict)
{
    _featureAngle = readScalar(snapDict.lookup("featureAngle"));
    _minEdgeForFeature = readLabel(snapDict.lookup("minEdgeForFeature"));
    _minFeatureEdgeLength = readScalar(snapDict.lookup("minFeatureEdgeLength"));
    _writeFeatures = readBool(snapDict.lookup("writeFeatures"));

    Info<< "  snapControls:"  << nl
        << "    - Feature angle              : " << _featureAngle  << nl
        << "    - Min edges for features     : " << _minEdgeForFeature << nl
        << "    - Min feature edge length    : " << _minFeatureEdgeLength << nl;

    const polyBoundaryMesh& bM = _polyMesh->boundaryMesh();
    const label NbPolyPatchs = bM.size();
    _triSurfList.resize(NbPolyPatchs, 0);
    _triSurfSearchList.resize(NbPolyPatchs, 0);
    _surfFeatList.resize(NbPolyPatchs, 0);
    _extEdgMeshList.resize(NbPolyPatchs, 0);
    _bndUseIntEdges.resize(NbPolyPatchs, true);
    _bndIsSnaped.resize(NbPolyPatchs, true);
    _bndLayers.resize(NbPolyPatchs);

    if (snapDict.found("boundaries"))
    {
        Info<< "    - Boundary specifications    : " << nl;

        const dictionary& bndDict = snapDict.subDict("boundaries");
        wordList bndDefined = bndDict.toc();

        forAll(bndDefined, patchI)
        {
            Info<< "    - " << bndDefined[patchI] << nl;

            const dictionary& patchDic = bndDict.subDict(bndDefined[patchI]);

            if (patchDic.found("triSurface"))
            {
                word file = patchDic.lookup("triSurface");

                bool exist = false;
                for(label patchJ = 0; patchJ < NbPolyPatchs && !exist; ++patchJ)
                {
                    exist = bM[patchJ].name() == bndDefined[patchI];

                    if (exist)
                    {
                        Info<< "        - Snaping surface    : " << file << nl;

                        _bndIsSnaped[patchJ] = false;

                        IOobject surfFile
                        (
                            file,
                            _polyMesh->time().constant(),
                            "triSurface",
                            _polyMesh->time(),
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        );

                        triSurface* bnd = new triSurface(surfFile.filePath());
                        addTriFace(patchJ, bnd);

                        if (patchDic.found("internalFeatureEdges"))
                        {
                            _bndUseIntEdges[patchJ] = readBool
                            (
                                patchDic.lookup("internalFeatureEdges")
                            );

                            Info<< "        - Use internal edges : "
                                << _bndUseIntEdges[patchJ] << nl;
                        }
                        else
                        {
                            _bndUseIntEdges[patchJ] = false;
                        }

                        if (patchDic.found("boundaryLayer"))
                        {
                            _bndLayers[patchJ] = SmootherBoundaryLayer
                            (
                                patchDic.subDict("boundaryLayer")
                            );
                        }
                    }
                }

                if (!exist)
                {
                    WarningIn("Foam::MeshSmoother::analyseDict()")
                        << "Patch " << bndDefined[patchI]
                        << " definied in smootherDict is not existing in "
                        << "polyMesh, existing patch in polyMesh ares: "
                        << bM.names() << nl;
                }
            }
        }
    }
    Info<< nl;
}

labelList Foam::SmootherBoundary::analyseFeatures
(
    List<labelHashSet> &pp,
    std::set<std::set<label> >& fP
)
{
    Info<< "  Analyse features" << nl;

    labelList pointType(_polyMesh->nPoints(), INTERIOR);

    forAll(_polyMesh->boundaryMesh(), patchI)
    {
        // points adressing
        std::map<label, label> p2s; // poly to tri surface
        std::map<label, label> s2p; // tri surface to poly

        const List<labelledTri> triFace = analyseBoundaryFace(patchI, p2s, s2p);

        geometricSurfacePatchList patchName;
        const polyPatch &patch = _polyMesh->boundaryMesh()[patchI];
        patchName.append(geometricSurfacePatch(word(""), patch.name(), patchI));

        // Renumber points
        pointField surfacePoints(p2s.size());
        std::map<label, label>::iterator ptI;
        for(ptI = p2s.begin(); ptI != p2s.end(); ++ptI)
        {
            surfacePoints[ptI->second] = _polyMesh->points()[ptI->first];
        }

        triSurface* triSurf = new triSurface(triFace, patchName, surfacePoints);
        surfaceFeatures* sF;
        sF = new surfaceFeatures(*triSurf, _featureAngle, 0, 0, false);

        sF->trimFeatures
        (
            _minFeatureEdgeLength,
            _minEdgeForFeature,
            _featureAngle
        );

        markPts(sF, s2p, _bndUseIntEdges[patchI], pointType, pp, fP);

        if (_triSurfSearchList[patchI] == 0)
        { // Store trisurface from blockMesh

            addTriFace(patchI, triSurf);
        }
    }

    // If use internal edges, mark feature points
    const scalar minCos = Foam::cos(degToRad(180.0 - _featureAngle));
    const pointField& polyPts = _polyMesh->points();
    forAll(pp, ptI)
    {
        if (pointType[ptI] == EDGE)
        {
            // Test if one of the feature surface don't use internal edges
            labelHashSet& featSet = _pointFeatureSet[ptI];
            bool canBeVertex = true;
            for
            (
                labelHashSet::iterator featI = featSet.begin();
                featI != featSet.end() && canBeVertex;
                ++featI
            )
            {
                canBeVertex = _bndUseIntEdges[featI.key()];
            }

            if (canBeVertex)
            {
                if (pp[ptI].size() == 2)
                { // Feature edge, check angle

                    DynamicList<vector> v(2);
                    forAllConstIter(labelHashSet, pp[ptI],neiI)
                    {
                        v.append(polyPts[neiI.key()] - polyPts[ptI]);
                        v.last() /= mag(v.last());
                    }

                    if (mag(v[0] & v[1]) < minCos)
                    {
                        pointType[ptI] = VERTEX;
                    }


                }
    //            Info<< ptI << " - ";
    //            forAllConstIter(labelHashSet, pp[ptI], ppI)
    //            {
    //                Info<< ppI.key() << " ";
    //            }
    //            Info<< nl;

                else if (pp[ptI].size() > 2)
                { // Feature point

                    pointType[ptI] = VERTEX;
                }
            }
        }
    }

    return pointType;
}

void Foam::SmootherBoundary::addTriFace
(
    const label patch,
    triSurface *triSurf
)
{
    _triSurfList[patch] = triSurf;
    _triSurfSearchList[patch] = new triSurfaceSearch(*triSurf);

    boolList surfBafReg(triSurf->patches().size());
    const polyBoundaryMesh& pBM = _polyMesh->boundaryMesh();
    forAll(triSurf->patches(), patchI)
    {
        surfBafReg[patchI] = (pBM[patchI].type() == "baffle");
    }

    surfaceFeatures* sF = new surfaceFeatures
    (
        *triSurf,
        _featureAngle,
        _minFeatureEdgeLength,
        _minEdgeForFeature,
        false
    );

    _extEdgMeshList[patch] = new extendedEdgeMesh(*sF, surfBafReg);
    _surfFeatList[patch] = sF;
}

Foam::List<Foam::labelledTri> Foam::SmootherBoundary::analyseBoundaryFace
(
    const label patchI,
    std::map<label, label> &p2s,
    std::map<label, label> &s2p
)
{
    List<labelledTri> triFaces;

    forAll(_polyMesh->boundaryMesh()[patchI], faceI)
    {
        const face& f = _polyMesh->boundaryMesh()[patchI][faceI];
        const label faceRef = getRefFromFace(f);

        labelList ptPoly(3);
        ptPoly[0] = _polyMesh->tetBasePtIs()[faceRef];
        ptPoly[1] = f.fcIndex(ptPoly[0]);

        for (label i = 2; i < f.size(); i++)
        {
            ptPoly[2] = f.fcIndex(ptPoly[1]);

            labelList ptTri(3);
            forAll(ptTri, ptI)
            {
                std::map<label, label>::iterator iter;
                iter = p2s.find(f[ptPoly[ptI]]);

                if (iter != p2s.end())
                {
                    ptTri[ptI] = iter->second;
                }
                else
                {
                    const label ref = p2s.size();
                    std::map<label,label>::iterator it = p2s.end();
                    p2s.insert(it , std::pair<label,label>(f[ptPoly[ptI]],ref));
                    s2p.insert(std::pair<label ,label>(ref, f[ptPoly[ptI]]));
                    ptTri[ptI] = ref;

                    _pointFeature[f[ptPoly[ptI]]] = patchI;
                    labelHashSet& featureSet = _pointFeatureSet[f[ptPoly[ptI]]];
                    featureSet.insert(patchI);
                }
            }

            triFaces.append(labelledTri(ptTri[0], ptTri[1], ptTri[2], patchI));
            ptPoly[1] = ptPoly[2];
        }
    }
    return triFaces;
}

Foam::label Foam::SmootherBoundary::getRefFromFace(const face &faceI)
{
    labelHashSet faceRef;
    faceRef.insert(_polyMesh->pointFaces()[faceI[0]]);
    forAll(faceI, ptI)
    {
        faceRef &= labelHashSet(_polyMesh->pointFaces()[faceI[ptI]]);
    }
    return faceRef.begin().key();
}

void Foam::SmootherBoundary::markPts
(
    surfaceFeatures* surfFeat,
    std::map<label,label> &s2p,
    const bool uE,
    labelList &pointType,
    List<labelHashSet>& pp,
    std::set<std::set<label> >& fP
)
{
    const triSurface& triSurf = surfFeat->surface();
    const edgeList& edgeLst = triSurf.edges();
    const labelList& featEdges = surfFeat->featureEdges();

    // Store faces points
    forAll(triSurf, faceI)
    {
        std::set<label> facePt;
        forAll(triSurf[faceI], ptI)
        {
            facePt.insert(s2p[triSurf[faceI][ptI]]);
        }
        fP.insert(facePt);
    }

    // Mark boundary points
    forAll(edgeLst, boundaryEdgI)
    {
        const edge& edg = edgeLst[boundaryEdgI];

        pointType[s2p[edg.start()]] = BOUNDARY;
        pointType[s2p[edg.end()]] = BOUNDARY;
    }

    // Create set of edges existing in polyMesh
    std::set<std::set<label> >  polyMeshEdges;
    forAll(_polyMesh->edges(), edgeI)
    {
        std::set<label> edgeDef;
        edgeDef.insert(_polyMesh->edges()[edgeI].first());
        edgeDef.insert(_polyMesh->edges()[edgeI].last());

        polyMeshEdges.insert(edgeDef);
    }

    // Mark feature edges points
    if (!uE)
    { // If use internal feature edge, mark boundary edges only

        forAll(featEdges, edgeI)
        {
            if (triSurf.edgeFaces()[featEdges[edgeI]].size() == 1)
            {
                const edge& edg = edgeLst[featEdges[edgeI]];

                const label pt1 = s2p[edg.start()];
                const label pt2 = s2p[edg.end()];

                // Test if edge is a polyMesh edge, dont know why but some time
                // the triSurface feature edge is not a polyMeshEdge..
                std::set<label> edgeDef;
                edgeDef.insert(pt1);
                edgeDef.insert(pt2);
                if (polyMeshEdges.find(edgeDef) != polyMeshEdges.end())
                {
                    pointType[pt1] = EDGE;
                    pointType[pt2] = EDGE;

                    pp[pt1].insert(pt2);
                    pp[pt2].insert(pt1);
                }
            }
        }
    }
    else
    {
        forAll(featEdges, edgeI)
        {
            const edge& edg = edgeLst[featEdges[edgeI]];

            const label pt1 = s2p[edg.start()];
            const label pt2 = s2p[edg.end()];

            std::set<label> edgeDef;
            edgeDef.insert(pt1);
            edgeDef.insert(pt2);
            if (polyMeshEdges.find(edgeDef) != polyMeshEdges.end())
            {
                pointType[pt1] = EDGE;
                pointType[pt2] = EDGE;

                pp[pt1].insert(pt2);
                pp[pt2].insert(pt1);
            }
        }
    }


//    // Mark vertex points TODO add edges points
//    const labelList& fPts = surfFeat->featurePoints();
//    forAll(fPts, vertexI)
//    {
//        pointType[s2p[fPts[vertexI]]] = VERTEX;
//    }
}

void Foam::SmootherBoundary::createPoints(labelList &pointType)
{
    label nbVertex = 0, nbEdge = 0, nbBoundary = 0, nbInterior = 0;

    forAll(pointType, ptI)
    {
        if (pointType[ptI] == VERTEX)
        {
            ++nbVertex;
            _point[ptI] = new SmootherVertex(ptI, _polyMesh->points()[ptI]);
            _featuresPoint.insert(ptI);
        }
        else if (pointType[ptI] == EDGE)
        {
            ++nbEdge;
            _point[ptI] = new SmootherEdge
            (
                ptI,
                _pointFeature[ptI],
                _polyMesh->points()[ptI]
            );

            if (!_bndIsSnaped[_pointFeature[ptI]])
            {
                _unsnapedPoint.insert(ptI);
            }
            _featuresPoint.insert(ptI);
        }
        else if (pointType[ptI] == BOUNDARY)
        {
            ++nbBoundary;
            _point[ptI] = new SmootherSurface
            (
                ptI,
                _pointFeature[ptI],
                _polyMesh->points()[ptI]
            );

            if (!_bndIsSnaped[_pointFeature[ptI]])
            {
                _unsnapedPoint.insert(ptI);
            }
            _featuresPoint.insert(ptI);
        }
        else if (pointType[ptI] == INTERIOR)
        {
            ++nbInterior;
            _point[ptI] = new SmootherPoint(ptI, _polyMesh->points()[ptI]);
            _interiorPoint.insert(ptI);
        }
    }

    Info<< "      - Number of feature points:  " << nbVertex << nl
        << "      - Number of edge points:     " << nbEdge << nl
        << "      - Number of boundary points: " << nbBoundary << nl
        << "      - Number of interior points: " << nbInterior << nl << nl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SmootherBoundary::SmootherBoundary
(
    dictionary &snapDict,
    polyMesh* mesh
)
:
    _polyMesh(mesh),
    _point(List<SmootherPoint*>(_polyMesh->nPoints()))
{
    analyseDict(snapDict);
    List<labelHashSet> pp(mesh->nPoints());
    std::set<std::set<label> > fP;
    labelList pointType = analyseFeatures(pp, fP);
    createPoints(pointType);

    if (_writeFeatures)
    {
        writeFeatures(pointType, pp, fP);
    }
}

// * * * * * * * * * * * * * * * * Desctructor  * * * * * * * * * * * * * * //

Foam::SmootherBoundary::~SmootherBoundary()
{
    forAll(_triSurfList, trisurfaceI)
    {
        delete _triSurfList[trisurfaceI];
    }

    forAll(_triSurfSearchList, triSurfSearchI)
    {
        delete _triSurfSearchList[triSurfSearchI];
    }

    forAll(_surfFeatList, surfFeatI)
    {
        delete _surfFeatList[surfFeatI];
    }

    forAll(_extEdgMeshList, extEdgMeshI)
    {
        delete _extEdgMeshList[extEdgMeshI];
    }

    forAll(_point, ptI)
    {
        delete _point[ptI];
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SmootherBoundary::writeFeatures
(
    labelList &pointType,
    List<labelHashSet> &pp,
    std::set<std::set<label> >& fP
) const
{
    const pointField& pt = _polyMesh->points();

    // ------------------------------------------------------------------------
    // Feature points
    labelHashSet vertexPt;
    std::ofstream vOutClean("featurePoints.vtk");
    vOutClean.close();
    std::ofstream vOut("featurePoints.vtk", std::ios::app);

    // Write header
    vOut<< "# vtk DataFile Version 2.0" << nl << "mesh vertex as vtk" << nl
        << "ASCII" << nl << nl << "DATASET POLYDATA" << nl;

    forAll(pointType, ptI)
    {
        if (pointType[ptI] == VERTEX)
        {
            vertexPt.insert(ptI);
        }
    }

    // Write points
    vOut<< "POINTS " << vertexPt.size() << " float" << nl;
    forAllConstIter(labelHashSet, vertexPt, ptI)
    {
        const label& ptR = ptI.key();
        vOut<< pt[ptR].x() << " " << pt[ptR].y() << " " << pt[ptR].z() << nl;
    }

    vOut.close();

    // ------------------------------------------------------------------------
    // Feature edges initial
    std::ofstream eOutClean("featureEdges.vtk");
    eOutClean.close();
    std::ofstream eOut("featureEdges.vtk", std::ios::app);

    // Map of points poly to vtk ref
    std::map<label, label> p2vtk;
    std::set<std::set<label> > edges;
    forAll(pointType, ptI)
    {
        if (pointType[ptI] == VERTEX || pointType[ptI] == EDGE)
        {
            p2vtk.insert(std::pair<label, label>(ptI, p2vtk.size()));

            forAllConstIter(labelHashSet, pp[ptI], ePtI)
            {
                std::set<label> edgePt;
                edgePt.insert(ptI);
                edgePt.insert(ePtI.key());

                edges.insert(edgePt);
            }
        }
    }

    // Write header
    eOut<< "# vtk DataFile Version 2.0" << nl << "mesh edges as vtk" << nl
        << "ASCII" << nl << nl << "DATASET POLYDATA" << nl;

    // Write points
    eOut<< "POINTS " << p2vtk.size() << " float" << nl;
    for
    (
        std::map<label, label>::iterator ptI = p2vtk.begin();
        ptI != p2vtk.end();
        ++ptI
    )
    {
        const label& ptR = ptI->first;
        eOut<< pt[ptR].x() << " " << pt[ptR].y() << " " << pt[ptR].z() << nl;
    }

    // Write lines
    eOut<< nl << "LINES " << edges.size() << " " << edges.size()*3 << nl;
    for
    (
        std::set<std::set<label> >::iterator edgeI = edges.begin();
        edgeI != edges.end();
        ++edgeI
    )
    {
        eOut<< 2 << " "
            << p2vtk.at(*(*edgeI).begin()) << " "
            << p2vtk.at(*(*edgeI).rbegin()) << nl;
    }

    eOut.close();

    // ------------------------------------------------------------------------
    // Boundary
    std::ofstream bOutClean("boundary.vtk");
    bOutClean.close();
    std::ofstream bOut("boundary.vtk", std::ios::app);

    std::map<label, label> p2vtk2;
    forAll(pointType, ptI)
    {
        if (pointType[ptI] == VERTEX || pointType[ptI] == EDGE || pointType[ptI] == BOUNDARY)
        {
            p2vtk2.insert(std::pair<label, label>(ptI, p2vtk2.size()));
        }
    }

    // Write header
    bOut<< "# vtk DataFile Version 2.0" << nl << "mesh boundaries as vtk" << nl
        << "ASCII" << nl << nl << "DATASET POLYDATA" << nl;

    // Write points
    bOut<< "POINTS " << p2vtk2.size() << " float" << nl;
    for
    (
        std::map<label, label>::iterator ptI = p2vtk2.begin();
        ptI != p2vtk2.end();
        ++ptI
    )
    {
        const label& ptR = ptI->first;
        bOut<< pt[ptR].x() << " " << pt[ptR].y() << " " << pt[ptR].z() << nl;
    }

    // Write triangles
    bOut<< nl << "POLYGONS " << fP.size() << " " << fP.size()*4 << nl;
    for
    (
        std::set<std::set<label> >::iterator faceI = fP.begin();
        faceI != fP.end();
        ++faceI
    )
    {
        bOut<< 3;

        for
        (
            std::set<label>::iterator ptI = faceI->begin();
            ptI != faceI->end();
            ++ptI
        )
        {
            bOut<< " " << p2vtk2.at(*ptI);
        }
        bOut<< nl;
    }

    bOut.close();

    // ------------------------------------------------------------------------
    // Interior points
    labelHashSet interiorPt;
    std::ofstream iOutClean("interiorPoints.vtk");
    iOutClean.close();
    std::ofstream iOut("interiorPoints.vtk", std::ios::app);

    // Write header
    iOut<< "# vtk DataFile Version 2.0" << nl << "mesh points as vtk" << nl
        << "ASCII" << nl << nl << "DATASET POLYDATA" << nl;

    forAll(pointType, ptI)
    {
        if (pointType[ptI] == INTERIOR)
        {
            interiorPt.insert(ptI);
        }
    }

    // Write points
    iOut<< "POINTS " << interiorPt.size() << " float" << nl;
    forAllConstIter(labelHashSet, interiorPt, ptI)
    {
        const label& ptR = ptI.key();
        iOut<< pt[ptR].x() << " " << pt[ptR].y() << " " << pt[ptR].z() << nl;
    }

    iOut.close();
}

void Foam::SmootherBoundary::removeSnapPoint(const label ref)
{
    _unsnapedPoint.erase(ref);
}

void SmootherBoundary::writeAllSurfaces(const label iterRef) const
{
    forAll(_triSurfSearchList, surfI)
    {
        std::ostringstream oss;
        oss << _triSurfSearchList[surfI]->surface().patches().begin()->name()
            << "-" << "Iter-" << iterRef
            << ".obj";
        _triSurfSearchList[surfI]->surface().write(oss.str());
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
