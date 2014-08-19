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

#ifndef BLOCKMESHSMOOTHER_H
#define BLOCKMESHSMOOTHER_H

#include "argList.H"

#include "blockMeshTopology.h"

#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class blockMeshSmoother Declaration
\*---------------------------------------------------------------------------*/

class blockMeshSmoother
{
    // Private data


        // Pointer to blockMesh
        blockMesh *blockMeshPtr_;

        // Dictionary
        dictionary dict_;

        // List of cell linked to points
        labelListList pointCells_;

        // list of points for each cell
        labelListList cellPoints_;

        // cell linked to points
        List<std::set<label> > pointCellsS_;

        // cell linked to cells
        List<std::set<label> > cellNeighbors_;

        // Sum of cell quality for each pt
        scalarList sumCellQuality_;

        // List of cell quality
        scalarList cellQuality_;

        // Topology of points
        blockMeshTopology pointTopology_;

        // Min mesh quality
        scalar minQuality_;

        // Mean mesh quality
        scalar meanQuality_;

        // Write intermediate meshs
        bool writeIntermediateMesh_;



    // Private member functions

        // Mean ratio
        void meshMeanRatio();

        // AddTransformedElementNodesAndWeights
        pointField addTransformedElementNodesAndWeights
        (
            scalarList &wj,
            std::set<label> &tp,
            const scalar &targetQual,
            const scalar &eTP
        );

        void addUntransformedElementNodesAndWeights
        (
            pointField &pi,
            scalarList &wj,
            std::set<label> &tn,
            const scalar &targetQual
        );

        void computeNewNodes
        (
            pointField &pi,
            scalarList &wj,
            std::set<label> &tn
        );

        pointField iterativeNodeRelaxation
        (
            pointField &pi,
            std::set<label> &tn,
            const scalarList &rT,
            labelList &relaxedCells
        );

        //- as copy (not implemented)
        blockMeshSmoother(const blockMeshSmoother&);

public:

    //- Runtime type information
//    TypeName("blockMeshSmoother");


    // Constructors

        //- Construct from blockMesh, dictionary and args
        blockMeshSmoother
        (
            blockMesh *block,
            dictionary &smootherDict,
            const argList &args
        );

    //- Destructor
        ~blockMeshSmoother();

    // Member functions

        //- Start smoothing
        void smoothing(const argList &args);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // BLOCKMESHSMOOTHER_H

// ************************************************************************* //
