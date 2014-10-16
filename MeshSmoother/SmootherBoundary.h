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

#include "SmootherBoundaryLayer.h"

#include <map>
#include <set>

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
            VERTEX,
            INTERIOR_BL1,
            INTERIOR_BL2,
            INTERIOR_BL3
        };

        // Pointer to polyMesh
        polyMesh* _polyMesh;

        // Searchable surfaces list
        List<triSurface*> _triSurfList;
        List<triSurfaceSearch*> _triSurfSearchList;

        // Searchable edge list
        List<surfaceFeatures*> _surfFeatList;
        List<extendedEdgeMesh*> _extEdgMeshList;

        // Parameter of patch
        boolList _bndUseIntEdges;
        boolList _bndIsSnaped;
        List<SmootherBoundaryLayer> _bndLayers;

        // Point and patch feature ref
        std::map<label, label> _pointFeature;
        std::map<label, labelHashSet> _pointFeatureSet;

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
        bool _writeFeatures;

    //- Private member functions

        void analyseDict(dictionary &snapDict);

        labelList analyseFeatures
        (
            List<labelHashSet>& pp,
            std::set<std::set<label> >& fP
        );

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
            const bool uE,
            labelList &pointType,
            List<labelHashSet>& pp,
            std::set<std::set<label> >& fP
        );

        void createPoints(labelList &pointType);

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

        // Get hash set of specific points
        const labelHashSet& unSnapedPoints() const {return _unsnapedPoint;}
        const labelHashSet& interiorPoints() const {return _interiorPoint;}
        const labelHashSet& featuresPoints() const {return _featuresPoint;}

        // Write edges as VTK points
        void writeFeatures
        (
            labelList& pointType,
            List<labelHashSet>& pp,
            std::set<std::set<label> >& fP
        ) const;

        void removeSnapPoint(const label ref);

        void writeAllSurfaces(const label iterRef) const;
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
