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

#include "cellSmoother.h"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::cellSmoother::cellSmoother(pointField &H)
    :
      points_(H)
{
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

Foam::scalar Foam::cellSmoother::tetrahedralMeanRatio
(
    const label &pt,
    const label &pt1,
    const label &pt2,
    const label &pt3
) const
{
    const Tensor<scalar> mA
    (
        points_[pt1] - points_[pt],
        points_[pt2] - points_[pt],
        points_[pt3] - points_[pt]
    );
    const scalar sigma(det(mA));
    if (sigma > VSMALL)
    {
        return 3.0*std::pow(sigma, 2.0/3.0)/magSqr(mA);
    }

    return 0.0;
}

Foam::pointField Foam::cellSmoother::dualOctahedron() const
{
    pointField oct(6);

    oct[0] = (points_[0] + points_[1] + points_[2] + points_[3])/4;
    oct[1] = (points_[0] + points_[4] + points_[5] + points_[1])/4;
    oct[2] = (points_[1] + points_[5] + points_[6] + points_[2])/4;
    oct[3] = (points_[2] + points_[6] + points_[7] + points_[3])/4;
    oct[4] = (points_[0] + points_[3] + points_[7] + points_[4])/4;
    oct[5] = (points_[4] + points_[7] + points_[6] + points_[5])/4;

    return oct;
}

Foam::pointField Foam::cellSmoother::dualOctahedronFaceCentroid
(
    const pointField &oct
) const
{
    pointField octC(8);
    octC[0] = (oct[0] + oct[1] + oct[4])/3.0;
    octC[1] = (oct[0] + oct[2] + oct[1])/3.0;
    octC[2] = (oct[0] + oct[3] + oct[2])/3.0;
    octC[3] = (oct[0] + oct[4] + oct[3])/3.0;
    octC[4] = (oct[5] + oct[4] + oct[1])/3.0;
    octC[5] = (oct[5] + oct[1] + oct[2])/3.0;
    octC[6] = (oct[5] + oct[2] + oct[3])/3.0;
    octC[7] = (oct[5] + oct[3] + oct[4])/3.0;

    return octC;
}

Foam::pointField Foam::cellSmoother::dualOctahedronNormals
(
    const pointField &oct
) const
{
    pointField octN(8);

    octN[0] = (oct[1] - oct[0]) ^ (oct[4] - oct[0]);
    octN[1] = (oct[2] - oct[0]) ^ (oct[1] - oct[0]);
    octN[2] = (oct[3] - oct[0]) ^ (oct[2] - oct[0]);
    octN[3] = (oct[4] - oct[0]) ^ (oct[3] - oct[0]);
    octN[4] = (oct[4] - oct[5]) ^ (oct[1] - oct[5]);
    octN[5] = (oct[1] - oct[5]) ^ (oct[2] - oct[5]);
    octN[6] = (oct[2] - oct[5]) ^ (oct[3] - oct[5]);
    octN[7] = (oct[3] - oct[5]) ^ (oct[4] - oct[5]);

    return octN;
}

Foam::pointField Foam::cellSmoother::tranformedHexahedron
(
    const scalar &cor,
    const Foam::pointField &octC,
    const Foam::pointField &octN
) const
{
    pointField Hp(8);
    forAll (Hp, ptI)
    {
        Hp[ptI] = octC[ptI] + cor/std::sqrt(mag(octN[ptI]))*octN[ptI];
    }

    return Hp;
}

Foam::point Foam::cellSmoother::centroidOfHex(const pointField &Hp) const
{
    point c(0.0, 0.0, 0.0);
    forAll (points_, ptI)
    {
        c += Hp[ptI];
    }
    return c/8.0;
}

Foam::scalar Foam::cellSmoother::ratioOfAvgLength
(
    Foam::pointField &Hp
) const
{
    return
    (
        (edgeAverageLength()/12.0)/(cellSmoother(Hp).edgeAverageLength()/12.0)
    );
}

Foam::scalar Foam::cellSmoother::edgeAverageLength() const
{
    return
    (
        mag(points_[0] - points_[1]) +
        mag(points_[1] - points_[2]) +
        mag(points_[2] - points_[3]) +
        mag(points_[3] - points_[0]) +
        mag(points_[0] - points_[4]) +
        mag(points_[1] - points_[5]) +
        mag(points_[2] - points_[6]) +
        mag(points_[3] - points_[7]) +
        mag(points_[4] - points_[5]) +
        mag(points_[5] - points_[6]) +
        mag(points_[6] - points_[7]) +
        mag(points_[7] - points_[4])
    )/12.0;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::cellSmoother::meanRatio() const
{
    scalar qn(tetrahedralMeanRatio(0, 3, 4, 1));
    if (qn < VSMALL)
    {
        return 0.0;
    }
    scalar qt(qn);
    qn = tetrahedralMeanRatio(1, 0, 5, 2);
    if (qn < VSMALL)
    {
        return 0.0;
    }
    qt += qn;
    qn = tetrahedralMeanRatio(2, 1, 6, 3);
    if (qn < VSMALL)
    {
        return 0.0;
    }
    qt += qn;
    qn = tetrahedralMeanRatio(3, 2, 7, 0);
    if (qn < VSMALL)
    {
        return 0.0;
    }
    qt += qn;
    qn = tetrahedralMeanRatio(4, 7, 5, 0);
    if (qn < VSMALL)
    {
        return 0.0;
    }
    qt += qn;
    qn = tetrahedralMeanRatio(5, 4, 6, 1);
    if (qn < VSMALL)
    {
        return 0.0;
    }
    qt += qn;
    qn = tetrahedralMeanRatio(6, 5, 7, 2);
    if (qn < VSMALL)
    {
        return 0.0;
    }
    qt += qn;
    qn = tetrahedralMeanRatio(7, 6, 4, 3);
    if (qn < VSMALL)
    {
        return 0.0;
    }
    qt += qn;

    return qt/8.0;
}

Foam::pointField Foam::cellSmoother::geometricTranform

(
    const scalar &cor
) const
{
    // Compute dual octahedron
    const pointField oct(dualOctahedron());

    // Compute new points
    pointField Hp(tranformedHexahedron
    (
        cor,
        dualOctahedronFaceCentroid(oct),
        dualOctahedronNormals(oct))
    );

    // Scaling of new points
    const scalar scalingfact(ratioOfAvgLength(Hp));

    const pointField C(8, centroidOfHex(Hp));

    // Return scaled new points
    return (C + scalingfact*(Hp - C));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
