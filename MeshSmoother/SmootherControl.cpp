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

#include "SmootherControl.h"

#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SmootherControl::SmootherControl(dictionary *smootherDict)
{
    // Get smoother parameters
    const dictionary& smoothDic = smootherDict->subDict("smoothControls");
    _maxIterations = readLabel(smoothDic.lookup("maxIterations"));
    _transformParam = readScalar(smoothDic.lookup("transformParameter"));
    _meanImprovTol = readScalar(smoothDic.lookup("meanImprovTol"));
    _maxMinCycleNoChange = readLabel(smoothDic.lookup("maxMinCycleNoChange"));
    _meanRelaxTable = readList<scalar>(smoothDic.lookup("meanRelaxationTable"));
    _minRelaxTable = readList<scalar>(smoothDic.lookup("minRelaxationTable"));
    _snapRelaxTable = readList<scalar>(smoothDic.lookup("snapRelaxationTable"));
    _ratioForMin = readScalar(smoothDic.lookup("ratioWorstQualityForMin"));

    if (*_meanRelaxTable.rbegin() > VSMALL)
    {
        _meanRelaxTable.append(0.0);
    }
    if (*_minRelaxTable.rbegin() > VSMALL)
    {
        _minRelaxTable.append(0.0);
    }

    if (*_snapRelaxTable.rbegin() > VSMALL)
    {
        _snapRelaxTable.append(0.0);
    }
//    if (*_snapRelaxTable.begin() < (1.0 - VSMALL))
//    {
//        scalarList relax(_snapRelaxTable.size() + 1);
//        relax[0] = 1.0;

//        forAll(_snapRelaxTable, relaxI)
//        {
//            relax[relaxI + 1] = _snapRelaxTable[relaxI];
//        }
//        _snapRelaxTable = relax;
//    }

    Info<< nl
        << "  smoothControls:"  << nl
        << "    - Max iterations             : " << _maxIterations  << nl
        << "    - Tranformation parameter    : " << _transformParam << nl
        << "    - Mean improvement tolerance : " << _meanImprovTol << nl
        << "    - Max ineffective iteration  : " << _maxMinCycleNoChange << nl
        << "    - Mean relaxation table      : " << _meanRelaxTable << nl
        << "    - Min relaxation table       : " << _minRelaxTable << nl
        << "    - Snap relaxation table      : " << _snapRelaxTable << nl
        << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
