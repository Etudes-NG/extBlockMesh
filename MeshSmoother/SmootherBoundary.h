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

#ifndef MESHSMOOTHERBOUNDARY_H
#define MESHSMOOTHERBOUNDARY_H

#include "List.H"
#include "boolList.H"
#include "labelledTri.H"
#include "point.H"
#include "HashSet.H"

#include "extendedEdgeMesh.H"
#include "triSurfaceMesh.H"
#include "searchableSurfaces.H"
#include "triSurfaceMesh.H"
#include "surfaceFeatures.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class polyMesh;
class SmootherPoint;

/*---------------------------------------------------------------------------*\
                    Class MeshSmootherBoundary Declaration
\*---------------------------------------------------------------------------*/

class SmootherBoundary
{
    //- Private data

        enum pointType
        {
            INTERIOR,
            BOUNDARY,
            EDGE,
            VERTEX
        };

        // Pointer to polyMesh
        polyMesh* _polyMesh;

        // Searchable surfaces list
        List<triSurface*> _triSurfList;
        List<triSurfaceSearch*> _triSurfSearchList;

        // Searchable edge list
        List<surfaceFeatures*> _surfFeatList;
        List<extendedEdgeMesh*> _extEdgMeshList;

        // Map of boundary name
        boolList _bndUseIntEdges;
        boolList _bndIsSnaped;

        // Point types and feature ref
        labelList _pointType;
        std::map<label, label> _pointFeature;

        // Point as SmootherPoints
        List<SmootherPoint*> _point;

        // Hash set of specific points
        labelHashSet _unsnapedPoint;
        labelHashSet _featuresPoint;
        labelHashSet _interiorPoint;

        // Inputs snapControls
        scalar _featureAngle;
        scalar _minFeatureEdgeLength;
        label _minEdgeForFeature;

    //- Private member functions

        void analyseDict(dictionary &snapDict);

        void analyseFeatures();

        void addTriFace(const label patch, triSurface *triSurf);

        List<labelledTri> analyseBoundaryFace
        (
            const label patchI,
            std::map<label, label>& p2s,
            std::map<label, label>& s2p
        );

        label getRefFromFace(const face& faceI);

        void markPts
        (
            surfaceFeatures *surfFeat,
            std::map<label, label> &s2p,
            const bool uE
        );

public:

    //- Constructors

        SmootherBoundary(dictionary& snapDict, polyMesh* mesh);

    //- Destructor
    ~SmootherBoundary();

    //- Member functions

        // Get snaped point
        inline point snapToSurf(const label r, const point &pt) const;
        inline point snapToEdge(const label eRef, const point &pt) const;

        SmootherPoint* pt(const label p) const {return _point[p];}

        void createPoints();

        // Get hash set of specific points
        const labelHashSet& unSnapedPoints() const {return _unsnapedPoint;}
        const labelHashSet& interiorPoints() const {return _interiorPoint;}
        const labelHashSet& featuresPoints() const {return _featuresPoint;}

        // Write edges as VTK points
        void writeFeaturesEdges() const;

        void removeSnapPoint(const label ref);
};

point SmootherBoundary::snapToSurf(const label r, const point &pt) const
{
    const indexedOctree<treeDataTriSurface>& t = _triSurfSearchList[r]->tree();
    const point snapPoint = t.findNearest(pt, 1e10).hitPoint();
//    const scalar dist = mag(pt - snapPoint);
    return snapPoint;
}

point SmootherBoundary::snapToEdge(const label eRef, const point &pt) const
{
    return _extEdgMeshList[eRef]->edgeTree().findNearest(pt, 1e10).hitPoint();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // MESHSMOOTHERBOUNDARY_H

// ************************************************************************* //
