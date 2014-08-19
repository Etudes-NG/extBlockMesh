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

#ifndef MESHSMOOTHERPOINT_H
#define MESHSMOOTHERPOINT_H

#include "MeshSmoother.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class MeshSmootherPoint Declaration
\*---------------------------------------------------------------------------*/

class SmootherPoint
{
    //- Private data

protected:

    //- Protected data

        // Pointers
        static polyMesh* _polyMesh;
        static SmootherBoundary* _bnd;
        static SmootherParameter* _param;

        // Point storage
        point _initialPt; // Point before iteration
        point _movedPt;   // Moved point without relaxation
        point _relaxedPt; // Relaxed point

        // Scalar stored in points
        scalar _averageQuality;
        scalar _weightingFactor;

        // Relaxation level
        label _relaxLevel;

        // PolyMesh point ref
        label _ptRef;

        bool _isTransformed;

public:
    //- Constructors

        //- Construct from point ref
        SmootherPoint(const label ref);

        // Empty constructor
        SmootherPoint();

    //- Destructor
    virtual ~SmootherPoint() {}

    //- Member functions

        // Set boundary and control ref
        void setStaticItems
        (
            SmootherBoundary *bnd,
            SmootherParameter *param,
            polyMesh* mesh
        );

        // Set/get and reset quality
        void setQuality(const scalar &qual) {_averageQuality = qual;}
        const scalar &avgQual() const {return _averageQuality;}

        // Reset point
        inline void GETMeReset();
        inline void laplaceReset();

        // GETMe
        inline void addWeight(const scalar &wei, const point& pt);
        inline void addWeight(const scalar &wei);
        const scalar &weightingFactor() const {return _weightingFactor;}

        // Get points
        const point& getMovedPoint() const {return _movedPt;}
        const point& getRelaxedPoint() const {return _relaxedPt;}
        const point& getInitialPoint() const {return _initialPt;}

        // Move point
        virtual void GETMeSmooth() {_movedPt /= _weightingFactor;}
        virtual void laplaceSmooth();
        virtual void snap() {}
        virtual void featLaplaceSmooth() {}

        // Get information about point
        virtual void needSnap(){}
        virtual bool isSurface() const {return false;}
        virtual bool isEdge() const {return false;}

        // Set and reset relaxation level
        void resetRelaxationLevel() {_relaxLevel = 0;}
        inline void addRelaxLevel(const scalarList& r);

        // Compute relaxed point
        inline void relaxPoint(const scalarList &r);
};

void SmootherPoint::GETMeReset()
{
    _isTransformed = false;
    _weightingFactor = 0.0;
    _movedPt = point(0.0, 0.0, 0.0);
    _initialPt = _relaxedPt;
}

void SmootherPoint::laplaceReset()
{
    _isTransformed = true;
    _initialPt = _relaxedPt;
}

void SmootherPoint::addWeight(const scalar &wei, const point &pt)
{
    _isTransformed = true;
    _weightingFactor += wei;
    _movedPt += pt*wei;
}

void SmootherPoint::addWeight(const scalar& wei)
{
    _isTransformed = true;
    _weightingFactor += wei;
    _movedPt += wei*_initialPt;
}

void SmootherPoint::addRelaxLevel(const scalarList &r)
{
    if ((_relaxLevel + 1) < r.size())
    {
        ++_relaxLevel;
    }
}

void SmootherPoint::relaxPoint(const scalarList &r)
{
    if (_isTransformed)
    {
        _relaxedPt =
            (1.0 - r[_relaxLevel])*_initialPt + r[_relaxLevel]*_movedPt;
    }
    else
    {
        _relaxedPt = _initialPt;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // MESHSMOOTHERPOINT_H

// ************************************************************************* //
