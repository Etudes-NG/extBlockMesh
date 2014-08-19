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

#include "dictionary.H"
#include "polyMesh.H"
#include "Time.H"
#include <OFstream.H>

#include "SmootherVertex.h"
#include "SmootherEdge.h"
#include "SmootherSurface.h"

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::SmootherBoundary::analyseDict(dictionary &snapDict)
{
    _featureAngle = readScalar(snapDict.lookup("featureAngle"));
    _minEdgeForFeature = readLabel(snapDict.lookup("minEdgeForFeature"));
    _minFeatureEdgeLength = readScalar(snapDict.lookup("minFeatureEdgeLength"));

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

    if (snapDict.found("boundarys"))
    {
        Info<< "    - Boundary specifications    : " << nl;

        const dictionary& bndDict = snapDict.subDict("boundarys");
        wordList bndDefined = bndDict.toc();

        forAll(bndDefined, patchI)
        {
            Info<< "    - " << bndDefined[patchI] << nl;

            const dictionary& patchDic = bndDict.subDict(bndDefined[patchI]);

            if (patchDic.found("triSurface"))
            {
                word file = patchDic.lookup("triSurface");
                Info<< "        - Snaping surface    : " << file << nl;

                IOobject surfFile
                (
                    file,
                    _polyMesh->time().constant(),
                    "triSurface",
                    _polyMesh->time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                bool exist = false;
                for(label patchJ = 0; patchJ < NbPolyPatchs && !exist; ++patchJ)
                {
                    exist = bM[patchJ].name() == bndDefined[patchI];

                    if (exist)
                    {
                        _bndIsSnaped[patchJ] = false;
                        triSurface* bnd = new triSurface(surfFile.filePath());
                        addTriFace(patchJ, bnd);
                        _bndUseIntEdges[patchJ] = readBool
                        (
                            patchDic.lookup("internalFeatureEdges")
                        );
                        Info<< "        - Use internal edges : "
                            << _bndUseIntEdges[patchJ] << nl;
                    }
                }

                if (!exist)
                {
                    WarningIn("Foam::MeshSmoother::analyseDict()")
                        << "Patch " << bndDefined[patchI]
                        << " definied in smootherDict is not existing in "
                        << "polyMesh, existing patch in polyMesh ares: "
                        << bM.names();
                }
            }
        }
    }
    Info<< nl;
}

void Foam::SmootherBoundary::analyseFeatures()
{
    Info<< "  Analyse features" << nl;

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

        markPts(sF, s2p, _bndUseIntEdges[patchI]);

        if (_triSurfSearchList[patchI] == 0)
        { // Store trisurface from blockMesh

            addTriFace(patchI, triSurf);
        }
    }
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

    surfaceFeatures* sF;
    sF = new surfaceFeatures(*triSurf, _featureAngle, 0, 0, false);
    sF->trimFeatures(_minFeatureEdgeLength, _minEdgeForFeature, _featureAngle);

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
    const bool uE
)
{
    const triSurface& triSurf = surfFeat->surface();
    const edgeList& edgeLst = triSurf.edges();
    const labelList& featEdges = surfFeat->featureEdges();

    // Mark boundary points
    forAll(edgeLst, boundaryEdgI)
    {
        const edge& edg = edgeLst[boundaryEdgI];

        _pointType[s2p[edg.start()]] = BOUNDARY;
        _pointType[s2p[edg.end()]] = BOUNDARY;
    }

    if (!uE)
    {
        forAll(featEdges, edgeI)
        {
            if (triSurf.edgeFaces()[featEdges[edgeI]].size() == 1)
            {
                const edge& edg = edgeLst[featEdges[edgeI]];

                _pointType[s2p[edg.start()]] = EDGE;
                _pointType[s2p[edg.end()]] = EDGE;
            }
        }
    }
    else
    {
        forAll(featEdges, edgeI)
        {
            const edge& edg = edgeLst[featEdges[edgeI]];

            _pointType[s2p[edg.start()]] = EDGE;
            _pointType[s2p[edg.end()]] = EDGE;
        }

        // Mark vertex points TODO add edges points
        const labelList& fPts = surfFeat->featurePoints();
        forAll(fPts, vertexI)
        {
            _pointType[s2p[fPts[vertexI]]] = VERTEX;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SmootherBoundary::SmootherBoundary
(
    dictionary &snapDict,
    polyMesh* mesh
)
:
    _polyMesh(mesh),
    _pointType(labelList(_polyMesh->nPoints(), INTERIOR)),
    _point(List<SmootherPoint*>(_polyMesh->nPoints()))
{
    analyseDict(snapDict);
    analyseFeatures();
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

void Foam::SmootherBoundary::createPoints()
{
    label nbVertex = 0, nbEdge = 0, nbBoundary = 0, nbInterior = 0;

    forAll(_pointType, ptI)
    {
        if (_pointType[ptI] == VERTEX)
        {
            ++nbVertex;
            _point[ptI] = new SmootherVertex(ptI);
            _featuresPoint.insert(ptI);
        }
        else if (_pointType[ptI] == EDGE)
        {
            ++nbEdge;
            _point[ptI] = new SmootherEdge(ptI, _pointFeature[ptI]);

            if (!_bndIsSnaped[_pointFeature[ptI]])
            {
                _unsnapedPoint.insert(ptI);
            }
            _featuresPoint.insert(ptI);
        }
        else if (_pointType[ptI] == BOUNDARY)
        {
            ++nbBoundary;
            _point[ptI] = new SmootherSurface(ptI, _pointFeature[ptI]);

            if (!_bndIsSnaped[_pointFeature[ptI]])
            {
                _unsnapedPoint.insert(ptI);
            }
            _featuresPoint.insert(ptI);
        }
        else if (_pointType[ptI] == INTERIOR)
        {
            ++nbInterior;
            _point[ptI] = new SmootherPoint(ptI);
            _interiorPoint.insert(ptI);
        }
    }
    writeFeaturesEdges();

    Info<< "      - Number of feature points:  " << nbVertex << nl
        << "      - Number of edge points:     " << nbEdge << nl
        << "      - Number of boundary points: " << nbBoundary << nl
        << "      - Number of interior points: " << nbInterior << nl << nl;
}
#include <sstream>
void SmootherBoundary::writeFeaturesEdges() const
{

    labelHashSet edgesPt;
    std::ofstream eOutClean("edgeDir.vtk");
    eOutClean.close();
    std::ofstream eOut("edgeDir.vtk", std::ios::app);

    // Write header
    eOut<< "# vtk DataFile Version 2.0" << nl
        << "blockMeshDict edges as vtk" << nl
        << "ASCII" << nl << nl
        << "DATASET POLYDATA" << nl;

    forAll(_pointType, ptI)
    {
        if (_pointType[ptI] == VERTEX || _pointType[ptI] == EDGE)
        {
            edgesPt.insert(ptI);
        }
    }

    eOut<< "POINTS " << edgesPt.size() << " float" << nl;
    forAllConstIter(labelHashSet, edgesPt, ptI)
    {
        eOut<< _polyMesh->points()[ptI.key()].x() << " "
            << _polyMesh->points()[ptI.key()].y() << " "
            << _polyMesh->points()[ptI.key()].z() << nl;
    }
    eOut.close();

//    forAll(_surfFeatList, featI)
//    {
//        std::ostringstream s;
//        s << "features" << featI;
//        _surfFeatList[featI]->writeObj(s.str());
//    }
}

void Foam::SmootherBoundary::removeSnapPoint(const label ref)
{
    _unsnapedPoint.erase(ref);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
