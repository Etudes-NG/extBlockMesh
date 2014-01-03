#ifndef CELLSMOOTHER_H
#define CELLSMOOTHER_H

#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class cellSmoother Declaration
\*---------------------------------------------------------------------------*/

class cellSmoother
{
    // Private data

        // Cell points
        pointField points_;

    // Private member functions

        //- tetrahedral mean ratio
        scalar tetrahedralMeanRatio
        (
            const label &pt,
            const label &pt1,
            const label &pt2,
            const label &pt3
        ) const;
public:
    // Constructors

        //- Construct from blockMesh and dictionary
        cellSmoother(const pointField &H);

    // Member functions

        //- meanRatio
        scalar meanRatio() const;

        // geometricTranform
        // FIXME decompose geometricTransfom in multiple methods (trans, scal,.)
        pointField geometricTranform(const scalar &cor) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // CELLSMOOTHER_H

// ************************************************************************* //
