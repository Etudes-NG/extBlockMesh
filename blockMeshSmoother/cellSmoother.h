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
        pointField &points_;

    // Private member functions

        //- tetrahedral mean ratio
        scalar tetrahedralMeanRatio
        (
            const label &pt,
            const label &pt1,
            const label &pt2,
            const label &pt3
        ) const;

        //- Hexahedron dual octahedron
        pointField dualOctahedron() const;

        //- Octahedron faces centroid
        pointField dualOctahedronFaceCentroid(const pointField &oct) const;

        //- Octahedron normal
        pointField dualOctahedronNormals(const pointField &oct) const;

        //- Transfomred hex
        pointField tranformedHexahedron
        (
            const scalar &cor,
            const pointField &octC,
            const pointField &octN
        ) const;

        //- Centroid of hex
        point centroidOfHex(const pointField &Hp) const;

        //- Ratio of length
        scalar ratioOfAvgLength(pointField &Hp) const;

        //- Edge average length
        scalar edgeAverageLength() const;

public:
    // Constructors

        //- Construct from blockMesh and dictionary
        cellSmoother(pointField &H);

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
