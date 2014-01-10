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

        // Set of mobiles points
        std::set<label> mobilPoints_;

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
        (pointField &pi,
            scalarList &wj);

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
            blockMesh &block,
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
