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

#ifndef MESHSMOOTHERCONTROL_H
#define MESHSMOOTHERCONTROL_H

#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class dictionary;

/*---------------------------------------------------------------------------*\
                        Class MeshSmootherControl Declaration
\*---------------------------------------------------------------------------*/

class SmootherControl
{

    //- Private data

        scalarList _meanRelaxTable;
        scalarList _minRelaxTable;
        scalarList _snapRelaxTable;
        scalar _transformParam;
        scalar _meanImprovTol;
        scalar _ratioForMin;
        label _maxMinCycleNoChange;
        label _maxIterations;

public:
    //- Constructors

        //- Construct from dictionary
        SmootherControl(dictionary *smootherDict);

    //- Member functions

        // Get iteration limit
        const label& maxIteration() const {return _maxIterations;}

        // Get improvement limit for mean cycle
        const scalar& meanImprovTol() const {return _meanImprovTol;}

        // Get max of cycles without change during min cycles
        const label& maxMinCycleNoChange() const {return _maxMinCycleNoChange;}

        // Get relaxation table
        const scalarList& minRelaxTable() const {return _minRelaxTable;}
        const scalarList& meanRelaxTable() const {return _meanRelaxTable;}
        const scalarList& snapRelaxTable() const {return _snapRelaxTable;}

        // Get transformation parameter
        const scalar &transformationParameter() const {return _transformParam;}

        const scalar& ratioForMin() const {return _ratioForMin;}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // MESHSMOOTHERCONTROL_H

// ************************************************************************* //
