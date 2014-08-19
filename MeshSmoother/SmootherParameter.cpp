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

#include "SmootherParameter.h"

#include "polyMesh.H"
#include "Time.H"

#include "MeshSmoother.h"

#include <cstdio>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SmootherParameter::SmootherParameter
(
    SmootherControl* control,
    polyMesh *poly
)
:
    _ctrl(control),
    _polyMesh(poly),
    _iterNb(0),
    _nbMovedPoints(0),
    _nbRelaxations(0),
    _updateTime(0.0),
    _totalTime(0.0),
    _transformTreshold(1.0),
    _noMinImproveCounter(0)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SmootherParameter::printHeaders() const
{
    Info<< "| Iteration | Mean qual | Min qual  |    Time   |Smooth type|"
            "Nb pts move| Nb relax  |Nb unsnaped|" << nl
        << "|-----------|-----------|-----------|-----------|-----------|"
            "-----------|-----------|-----------|"<< nl;
}

void Foam::SmootherParameter::printStatus(const label nbUnsnaped)
{
    _updateTime = _polyMesh->time().elapsedCpuTime() - _updateTime;
    _totalTime += _updateTime;

    const char* cycle =
        (_actualCycle == meanCycleRunning)
        ?
            "mean"
            :
            (_actualCycle == snapCycleRunning)
            ?
                "snap" :
                "min ";

    std::printf
    (
        "|    %6i |   %6.4f  |   %6.4f  |  %6.2f   |    %s   |  %8i |    %6i "
        "|    %6i |\n",
        _iterNb,
        _meanQuality,
        _minQuality,
        _updateTime,
        cycle,
        _nbMovedPoints,
        _nbRelaxations,
        nbUnsnaped
    );
}

void Foam::SmootherParameter::printStats() const
{
    Info<< "============================================================="
           "====================================" << nl;
    const bool isConv = _iterNb < _ctrl->maxIteration() + 1;
    const char* conv = (isConv) ? "Converged in " : "Not converged";
    std::printf("%s %.3f s\n", conv,_totalTime);
    Info<< "=====================" << nl;
}

bool Foam::SmootherParameter::setSmoothCycle
(
    const scalar &meanImp,
    const scalar &minImp,
    const bool asUnSnaped,
    MeshSmoother* meshSmoother
)
{
    if (_actualCycle == snapCycleRunning)
    {
        if (asUnSnaped)
        {
            _actualCycle = meanCycleRunning;
        }

        if (_iterNb == _ctrl->maxIteration())
        {
            return false;
        }
        ++_iterNb;
        return true;
    }

    const scalar meanImprove = _meanQuality - meanImp;
    const scalar minImprove = _minQuality - minImp;
    const scalar& improvTol = _ctrl->meanImprovTol();

    _prevCycle = _actualCycle;
    if (_actualCycle == meanCycleRunning && meanImprove < improvTol)
    {
        _actualCycle = minCycleStart;
    }

    if (_actualCycle == minCycleRunning)
    {
        if (minImprove < VSMALL)
        { // No improvement
            ++_noMinImproveCounter;
        }

        if(_noMinImproveCounter > _ctrl->maxMinCycleNoChange())
        { // More than x time with no evolution
            _actualCycle = minCycleStart;
        }
    }

    if (_actualCycle == minCycleStart)
    {
        if(minImprove < improvTol && _prevCycle == minCycleRunning)
        {
            return false;
        }
        _actualCycle = minCycleRunning;
        _noMinImproveCounter = 0;

//        _transformTreshold = _minQuality + (1.0 - _minQuality)/2.0;
        _transformTreshold = meshSmoother->getTransformationTreshold();
    }

    if (_iterNb == _ctrl->maxIteration())
    {
        return false;
    }

    ++_iterNb;
    return true;
}

void Foam::SmootherParameter::setSmoothCycle(const bool asUnSnaped)
{
    if (asUnSnaped)
    {
        _actualCycle = snapCycleRunning;
    }
    else
    {
        _actualCycle = meanCycleRunning;
    }
    _prevCycle = _actualCycle;
}



void Foam::SmootherParameter::resetUpdateTime()
{
    _updateTime = _polyMesh->time().elapsedCpuTime();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
