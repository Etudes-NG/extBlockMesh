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

#ifndef MESHSMOOTHER_H
#define MESHSMOOTHER_H

#include "fvCFD.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class SmootherCell;
class polyMesh;
class blockMesh;
class SmootherControl;
class SmootherParameter;
class SmootherBoundary;

/*---------------------------------------------------------------------------*\
                      Class blockMeshSmoother Declaration
\*---------------------------------------------------------------------------*/

class MeshSmoother
{
    //- Private data

        // Pointer of parents
        polyMesh *_polyMesh;
        blockMesh *_blocks;

        // MeshSmoother child items
        SmootherControl* _ctrl;
        SmootherParameter* _param;
        SmootherBoundary* _bnd;

        // Smoother cell and points
        List<SmootherCell*> _cell;

    //- Private member functions

        // Quality analysis
        void analyseMeshQuality();
        void analyseMeshQuality(const labelHashSet &cell);
        void qualityStats();

        // GETMe smoothing
        labelHashSet addTransformedElementNodeWeight();
        void addUnTransformedElementNodeWeight(labelHashSet &tp);
        bool untransformedAndhavePointTransformed
        (
            const label cellI,
            const labelHashSet& tp
        );

        void iterativeNodeRelaxation(labelHashSet &tP, const scalarList &r);

        // Run one iteration
        bool runIteration();
        pointField getMovedPoints() const;

        // Write the mesh and meshQuality
        void writeMesh(const fvMesh& meshFv, volScalarField &meshQuality) const;

        // Smoothing algo
        void GETMeSmoothing();
        void snapSmoothing();

public:

    //- Constructors

        //- Construct from polyMesh and dictionary
        MeshSmoother
        (
            polyMesh *mesh,
            dictionary *smootherDict,
            blockMesh *blocks = 0
        );

    //- Destructor
    ~MeshSmoother();

    //- Member functions

        // Smooth the mesh acording to smoothDict
        void update();

        // Smooth the mesh acording to smoothDict and write intermediate mesh
        void updateAndWrite
        (
            word &regionName,
            word &defaultFacesName,
            word &defaultFacesType,
            Time &runTime
        );

        // Get tranformation treshold
        scalar getTransformationTreshold() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // MESHSMOOTHER_H

// ************************************************************************* //
