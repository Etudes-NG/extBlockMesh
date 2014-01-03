#include "cellSmoother.h"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellSmoother::cellSmoother(const pointField &H)
    :
      points_(H)
{
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::cellSmoother::meanRatio() const
{
    // Labels for quality computation
    labelList v1(8), v2(8), v3(8);
    v1[0] = 3; v1[1] = 0; v1[2] = 1; v1[3] = 2;
    v1[4] = 7; v1[5] = 4; v1[6] = 5; v1[7] = 6;

    v2[0] = 4; v2[1] = 5; v2[2] = 6; v2[3] = 7;
    v2[4] = 5; v2[5] = 6; v2[6] = 7; v2[7] = 4;

    v3[0] = 1; v3[1] = 2; v3[2] = 3; v3[3] = 0;
    v3[4] = 0; v3[5] = 1; v3[6] = 2; v3[7] = 3;

    scalar cQa(0.0);
    forAll (points_, ptI)
    {
        const point p1(points_[v1[ptI]] - points_[ptI]);
        const point p2(points_[v2[ptI]] - points_[ptI]);
        const point p3(points_[v3[ptI]] - points_[ptI]);
        const Tensor<scalar> mA(p1, p2, p3);

        const scalar sigma(det(mA));

        if (sigma > VSMALL)
        {
            cQa += 3*std::pow(sigma, 2.0/3.0)/magSqr(mA);
        }
    }

    return cQa/8.0;
}

Foam::pointField Foam::cellSmoother::geometricTranform
(
    const scalar &cor
) const
{
    // Labels for dual octahedron
    labelList vb1(6),  vb2(6), vb3(6), vb4(6);
    vb1[0] = 0;	vb1[1] = 0;	vb1[2] = 1;	vb1[3] = 2;	vb1[4] = 0;	vb1[5] = 4;
    vb2[0] = 1;	vb2[1] = 4;	vb2[2] = 5;	vb2[3] = 6;	vb2[4] = 3;	vb2[5] = 7;
    vb3[0] = 2;	vb3[1] = 5;	vb3[2] = 6;	vb3[3] = 7;	vb3[4] = 7;	vb3[5] = 6;
    vb4[0] = 3;	vb4[1] = 1;	vb4[2] = 2;	vb4[3] = 3;	vb4[4] = 4;	vb4[5] = 5;

    // Labels for normals
    labelList vc1(8),  vc2(8), vc3(8);
    vc1[0] = 0;	vc1[1] = 0;	vc1[2] = 0;	vc1[3] = 0;
    vc1[4] = 5;	vc1[5] = 5;	vc1[6] = 5;	vc1[7] = 5;

    vc2[0] = 1;	vc2[1] = 2;	vc2[2] = 3;	vc2[3] = 4;
    vc2[4] = 4;	vc2[5] = 1;	vc2[6] = 2;	vc2[7] = 3;

    vc3[0] = 4;	vc3[1] = 1;	vc3[2] = 2;	vc3[3] = 3;
    vc3[4] = 1;	vc3[5] = 2;	vc3[6] = 3;	vc3[7] = 4;

    // Compute dual octahedron
    pointField oct(6);
    forAll (vb1, octPtI)
    {
        oct[octPtI] =
        (
            points_[vb1[octPtI]] + points_[vb2[octPtI]] +
            points_[vb3[octPtI]] + points_[vb4[octPtI]]
        )/4;
    }

    // Compute centroid of octahedron faces
    pointField octC(8);
    forAll (octC, octCI)
    {
        octC[octCI] = (oct[vc1[octCI]] + oct[vc2[octCI]] + oct[vc3[octCI]])/3;
    }

    // Compute normal of octahedron faces
    pointField octN(8);
    forAll (octN, ptI)
    {
        octN[ptI] =
                (oct[vc2[ptI]] - oct[vc1[ptI]]) ^
                (oct[vc3[ptI]] - oct[vc1[ptI]]);
    }

    // Compute new points
    pointField Hp(8);
    forAll (points_, ptI)
    {
        Hp[ptI] = octC[ptI] + cor/std::sqrt(mag(octN[ptI]))*octN[ptI];
    }

    // Scaling of new points
    // Centroid of cell
    point c(0, 0, 0);
    forAll (points_, ptI)
    {
        c += Hp[ptI];
    }
    c /= 8;

    // Labels for edge lengh
    labelList vd1(12),  vd2(12);
    vd1[0] = 0;	vd1[1] = 1;	vd1[2] = 2;	vd1[3] = 3;
    vd1[4] = 0;	vd1[5] = 1;	vd1[6] = 2;	vd1[7] = 3;
    vd1[8] = 4;	vd1[9] = 5;	vd1[10]= 6; vd1[11]= 7;

    vd2[0] = 1;	vd2[1] = 2;	vd2[2] = 3;	vd2[3] = 0;
    vd2[4] = 4;	vd2[5] = 5;	vd2[6] = 6;	vd2[7] = 7;
    vd2[8] = 5;	vd2[9] = 6;	vd2[10]= 7;	vd2[11]= 4;

    // Scaling factor (keeping avg edge lenght)
    scalar mh1(0.0), mh2(0.0);
    forAll (vd1, edgeI)
    {
        mh1 += mag(points_[vd1[edgeI]] - points_[vd2[edgeI]]);
        mh2 += mag(Hp[vd1[edgeI]] - Hp[vd2[edgeI]]);
    }
    mh1 /= 12;
    mh2 /= 12;
    const scalar scalingfact(mh1/mh2);

    pointField C(8, c);

    // Return scaled new points
    return C + scalingfact*(Hp - C);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
