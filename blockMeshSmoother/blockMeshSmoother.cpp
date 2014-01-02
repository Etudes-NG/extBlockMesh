#include "blockMeshSmoother.h"

#include "cellSmoother.h"

#include <algorithm>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::blockMeshSmoother::blockMeshSmoother
(
    blockMesh &block,
    dictionary &smootherDict
)
    :
      blockMeshPtr_(&block),
      dict_(smootherDict),
      pointCells_(labelListList(blockMeshPtr_->points().size())),
      cellPoints_(labelListList(blockMeshPtr_->cells().size())),
      pointCellsS_(List<std::set<label> >(blockMeshPtr_->points().size())),
      cellNeighbors_(List<std::set<label> >(blockMeshPtr_->cells().size())),
      sumCellQuality_(scalarList(blockMeshPtr_->points().size())),
      cellQuality_(scalarList(blockMeshPtr_->cells().size()))
{
    forAll (blockMeshPtr_->cells(), cellI)
    {
        cellPoints_[cellI] = blockMeshPtr_->cells()[cellI].pointsLabel();

        // Store points conectivity
        forAll (cellPoints_[cellI], ptI)
        {
            pointCells_[cellPoints_[cellI][ptI]].append(cellI);
            pointCellsS_[cellPoints_[cellI][ptI]].insert(cellI);
        }
    }
    forAll (blockMeshPtr_->points(), ptI)
    {
        for
        (
            std::set<label>::iterator iter = pointCellsS_[ptI].begin();
            iter != pointCellsS_[ptI].end();
            ++iter
        )
        {
            cellNeighbors_[*iter].insert(ptI);
        }
    }

    std::set<label> allPoints;
    label val(0);
    forAll (blockMeshPtr_->points(), pointI)
    {
        allPoints.insert(allPoints.end(),val);
        ++val;
    }

    // Create a set of fixed points
    std::set<label> fixedPoints;
    forAll (blockMeshPtr_->patches(), patchI)
    {
        forAll (blockMeshPtr_->patches()[patchI], faceI)
        {
            const labelList facePoints(
                        blockMeshPtr_->patches()[patchI][faceI].pointsLabel());
            forAll (facePoints, pointI)
            {
                fixedPoints.insert(facePoints[pointI]);
            }
        }
    }

    std::set_difference
    (
        allPoints.begin(),
        allPoints.end(),
        fixedPoints.begin(),
        fixedPoints.end(),
        std::inserter(mobilPoints_, mobilPoints_.begin())
                );
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::blockMeshSmoother::meshMeanRatio
(
    scalar &min,
    scalar &avg
)
{
    // Reset sumCellQuality_
    sumCellQuality_ = scalarList(blockMeshPtr_->points().size(), 0.0);
    min = 1.0;
    avg = 0;
    forAll (blockMeshPtr_->cells(), cellI)
    {
        cellQuality_[cellI] =
        cellSmoother
        (
            blockMeshPtr_->cells()[cellI].points
            (
             blockMeshPtr_->points()
            )
        ).meanRatio();

        if (cellQuality_[cellI] < min)
        {
            min = cellQuality_[cellI]; // Store the new qMin
        }
        avg += cellQuality_[cellI];

        forAll(cellPoints_[cellI], pointI)
        {
            sumCellQuality_[cellPoints_[cellI][pointI]] += cellQuality_[cellI];
        }
    }
    avg /= blockMeshPtr_->cells().size();
}

Foam::pointField Foam::blockMeshSmoother::addTransformedElementNodesAndWeights
(
        scalarList &wj,
        std::set<label> &tp,
        const scalar &targetQual,
        const scalar &eTP
)
{
    pointField pi(blockMeshPtr_->points().size(), point(0.0, 0.0, 0.0));

    forAll (blockMeshPtr_->cells(), cellI)
    {
        if (cellQuality_[cellI] < targetQual)
        {
            const pointField H
            (
                blockMeshPtr_->cells()[cellI].points
                (
                    blockMeshPtr_->points()
                )
            );

            const pointField Hr(cellSmoother(H).geometricTranform(eTP));

            // Add Transformed Element Nodes And Weights
            forAll (cellPoints_[cellI], pointI)
            {
                // compute the associated weight
                const scalar wji
                (
                    std::sqrt(sumCellQuality_[cellPoints_[cellI][pointI]] /
                    (cellNeighbors_[cellI].size()*cellQuality_[cellI]))
                );

                // add the weight to temporary weighted sum
                wj[cellPoints_[cellI][pointI]] += wji;

                // add the resul weighted nodes to temporary weighted pt
                pi[cellPoints_[cellI][pointI]] += Hr[pointI] * wji;

                // add point to set of tranformed points
                tp.insert(cellPoints_[cellI][pointI]);
            }
        }
    }

    return pi;
}

std::set<Foam::label> Foam::blockMeshSmoother::getTransformedAndFreePoints
(
    std::set<label> &tp
)
{
    std::set<label> tn; // Transformed Nodes
    std::set_intersection
    (
        tp.begin(),
        tp.end(),
        mobilPoints_.begin(),
        mobilPoints_.end(),
        std::inserter(tn, tn.begin())
    );
    return tn;
}

void Foam::blockMeshSmoother::addUntransformedElementNodesAndWeights
(
        pointField &pi,
        scalarList &wj,
        std::set<label> &tn,
        const scalar &targetQual
)
{
    forAll (blockMeshPtr_->cells(), cellI)
    {
        if (cellQuality_[cellI] > targetQual)
        { // Unstransformed cell
            forAll (cellPoints_[cellI], pointI)
            {
                if (tn.find(cellPoints_[cellI][pointI]) != tn.end())
                { // Unstranformed & have point transformed
                    // compute the associated weight
                    const label pt(cellPoints_[cellI][pointI]);
                    const scalar wji
                    (
                        std::sqrt(sumCellQuality_[pt] /
                        (cellNeighbors_[cellI].size()*cellQuality_[cellI]))
                    );

                    // add the weight to temporary weighted sum
                    wj[pt] += wji;

                    // add the weighted nodes to temporary weighted pt
                    pi[pt] += blockMeshPtr_->points()[pt] * wji;
                }
            }
        }
    }
}

void Foam::blockMeshSmoother::computeNewNodes
(
    pointField &pi,
    scalarList &wj,
    std::set<label> &tn
)
{
    for
    (
        std::set<label>::iterator ptI = tn.begin();
        ptI != tn.end();
        ++ptI
    )
    {
        if (wj[*ptI] > VSMALL)
        {
            pi[*ptI] /= wj[*ptI];
        }
        else
        {
            pi[*ptI] = blockMeshPtr_->points()[*ptI];
        }
    }
}

Foam::pointField Foam::blockMeshSmoother::iterativeNodeRelaxation
(
    pointField &pi,
    std::set<label> &tn,
    const scalarList &rT
)
{
    pointField pip(blockMeshPtr_->points()); // Storage of relaxed nodes
    labelList nR(blockMeshPtr_->points().size(), 0);
    label totalReversedCells(0), nbRelaxations(0);
    while (!tn.empty())
    {
        for
        (
            std::set<label>::iterator ptI = tn.begin();
            ptI != tn.end();
            ++ptI
        )
        {
            const scalar r(rT[nR[*ptI]]);
            pip[*ptI] = (1.0 - r)*blockMeshPtr_->points()[*ptI] + r*pi[*ptI];
        }

        // Test new cells
        tn.clear();
        forAll (blockMeshPtr_->cells(), cellI)
        {
            pointField Ht(cellPoints_[cellI].size());
            forAll (cellPoints_[cellI], ptI)
            {
                Ht[ptI] = pip[cellPoints_[cellI][ptI]];
            }

            if (cellSmoother(Ht).meanRatio() < VSMALL)
            { // Cell is invalid
                forAll (cellPoints_[cellI], ptI)
                { // Insert cell point in point invalid list
                    tn.insert(cellPoints_[cellI][ptI]);
                }
            }
        }
        std::set<label> tn2; // Transformed Nodes
        std::set_intersection
        (
            tn.begin(),
            tn.end(),
            mobilPoints_.begin(),
            mobilPoints_.end(),
            std::inserter(tn2, tn2.begin())
        );
        tn = tn2;

        label reversedCells(0);
        for
        (
            std::set<label>::iterator ptI = tn.begin();
            ptI != tn.end();
            ++ptI
        )
        { // For each invalid node
            if(nR[*ptI] < (rT.size() - 1))
            {
                ++nR[*ptI];
            }
            ++reversedCells;
        }
        totalReversedCells += reversedCells;
        ++nbRelaxations;
    }

    if (totalReversedCells != 0)
    {
        Info<< "    Reversed: " << totalReversedCells << " cells in "
            << nbRelaxations << " relaxations\n";
    }

    return pip;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockMeshSmoother::smoothing()
{
    const scalar eTP(readScalar(dict_.lookup("elemTransformParameter")));

    scalar qualityAvg, qualityMin;
    meshMeanRatio(qualityMin, qualityAvg);

    Info<< "\nStart GETMe smoothing with:"  << nl
        << "  - Max itrerations           : "
        << readLabel(dict_.lookup("maxGETMeIter")) << nl
        << "  - Tranformation parameter   : "
        << eTP << nl
        << "  - Improvement tolerance     : "
        << readScalar(dict_.lookup("improvementTolerance")) << nl
        << "  - Max ineffective iteration : "
        << readScalar(dict_.lookup("maxIneffectiveIteration")) << nl
        << "  - Mean relaxation table     : "
        << readList<scalar>(dict_.lookup("qMeanRelaxationTable")) << nl
        << "  - Min relaxation table      : "
        << readList<scalar>(dict_.lookup("qMinRelaxationTable")) << nl
        << "\nAnalysis before start: avg quality = " << qualityAvg
        << ", min quality = " << qualityMin << endl;

    label nbIterations(0), noMinImproveCounter(0), cycle(0);

    scalar targetQual(1.0), deltaMeanImprov(1.0), deltaMinImprov(1.0);
    scalar previousMean(qualityAvg), previousMin(qualityMin);

    scalarList rT = readList<scalar>(dict_.lookup("qMeanRelaxationTable"));

    bool converged(false);
    while
    (
        nbIterations < readLabel(dict_.lookup("maxGETMeIter")) &&
        !converged
    )
    {
        // Storage of Temporary Nodes And Weights sum
        scalarList wj(blockMeshPtr_->points().size(), 0.0);
        std::set<label> tp; // List of transformed nodes

        pointField pi
        (
            addTransformedElementNodesAndWeights
            (
                wj,
                tp,
                targetQual,
                eTP
            )
        );

        std::set<label> tn(getTransformedAndFreePoints(tp));

        addUntransformedElementNodesAndWeights(pi, wj, tn, targetQual);

        computeNewNodes(pi, wj, tn);

        const pointField pip(iterativeNodeRelaxation(pi, tn, rT));

        // Update the mesh with new points
        // FIXME change setPoint to change all points
        forAll (blockMeshPtr_->points(), ptI)
        {
            blockMeshPtr_->setPoint(ptI, pip[ptI]);
        }

        // Compute new min and avg quality
        meshMeanRatio(qualityMin, qualityAvg);

        deltaMeanImprov = qualityAvg - previousMean;
        previousMean = qualityAvg;

        deltaMinImprov = qualityMin - previousMin;
        previousMin = qualityMin;

        if
        (
            cycle == 0 &&
            deltaMeanImprov < readScalar(dict_.lookup("improvementTolerance"))
        )
        {
            cycle = 1;
            rT = readList<scalar>(dict_.lookup("qMinRelaxationTable"));
        }

        if (cycle == 0)
        {
            Info<< "mean cycle ";
        }

        if (cycle == 1)
        {
            if (previousMin >= qualityMin)
            { // No improvement
                ++noMinImproveCounter;
            }
            if
            (
                noMinImproveCounter > readScalar
                (
                    dict_.lookup("maxIneffectiveIteration")
                )
            )
            {
                cycle = 2;
            }
            else
            {
                Info<< "min cycle  ";
            }
        }

        if (cycle == 2)
        {
            if
            (
                deltaMinImprov < readScalar
                (
                    dict_.lookup("improvementTolerance")
                )
            )
            {
                converged = true;
            }
            cycle = 1;
            noMinImproveCounter = 0;

            scalarList cqs(cellQuality_);
            std::sort(cqs.begin(), cqs.end());

            targetQual = cqs[blockMeshPtr_->cells().size()*0.5];

            Info<< "converged  ";
        }

        Info<< "  iter " << nbIterations << "\t avg quality = "
            << qualityAvg << "\t min quality = " << qualityMin << endl;

        ++nbIterations;
    }
    if (nbIterations != readLabel(dict_.lookup("maxGETMeIter")))
    {
        Info<< "GETMe smoothing converged in " << nbIterations
            << " iterations" << endl;
    }
    else
    {
        Info<< "GETMe smoothing not converged!" << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
