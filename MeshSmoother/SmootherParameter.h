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

#ifndef MESHSMOOTHERPARAMETER_H
#define MESHSMOOTHERPARAMETER_H

#include "scalarList.H"

#include "SmootherControl.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class polyMesh;
class MeshSmoother;

/*---------------------------------------------------------------------------*\
                   Class MeshSmootherParameter Declaration
\*---------------------------------------------------------------------------*/

class SmootherParameter
{
    //- Private data

        // Enum iteration type
        enum cycleStatus
        {
            snapCycleRunning,
            selectedSnapCycleRunning,
            meanCycleRunning,
            minCycleStart,
            minCycleRunning
        };

        // Pointer of parents
        SmootherControl* _ctrl;
        polyMesh *_polyMesh;

        // Parameters
        label _iterNb;
        label _actualCycle;
        label _prevCycle;
        label _nbMovedPoints;
        label _nbRelaxations;
        scalar _updateTime;
        scalar _totalTime;
        scalar _transformTreshold;
        label _noMinImproveCounter;
        scalar _minQuality;
        scalar _meanQuality;

public:
    //- Constructors

        //- Construct from polyMesh and dictionary
        SmootherParameter(SmootherControl *control, polyMesh* poly);

    //- Member functions

        // Print
        void printHeaders() const;
        void printStatus(const label nbUnsnaped);
        void printStats() const;

        // Change smooth cycle
        bool setSmoothCycle
        (
            const scalar &meanImp,
            const scalar &minImp,
            const bool asUnSnaped,
            MeshSmoother* meshSmoother
        );
        void setSmoothCycle(const bool asUnSnaped);

        // Get iter number
        const label &getIterNb() const {return _iterNb;}
        void setIterNb() {_iterNb = 1;}

        // Set/get quality stats
        void setMinQual(const scalar& qual) {_minQuality = qual;}
        const scalar& minQual() const {return _minQuality;}

        void setMeanQual(const scalar& qual) {_meanQuality = qual;}
        const scalar &meanQual() const {return _meanQuality;}

        const scalar &transformationTreshold() const{return _transformTreshold;}
        inline const scalarList& relaxationTable() const;

        void setNbMovedPoints(const scalar& nbMoved) {_nbMovedPoints = nbMoved;}
        void setNbRelaxations(const label nbRelax) {_nbRelaxations = nbRelax;}

        void resetUpdateTime();
};

const scalarList &SmootherParameter::relaxationTable() const
{
    if (_actualCycle == meanCycleRunning)
    {
        return _ctrl->meanRelaxTable();
    }
    else
    {
        return _ctrl->minRelaxTable();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // MESHSMOOTHERPARAMETER_H

// ************************************************************************* //
