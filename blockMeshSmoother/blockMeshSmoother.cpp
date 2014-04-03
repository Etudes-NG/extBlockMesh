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

#include "blockMeshSmoother.h"

#include "Time.H"
#include "cellSmoother.h"

#include <algorithm>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::blockMeshSmoother::blockMeshSmoother
(
    blockMesh *block,
    dictionary &smootherDict,
    const argList &args
)
:
    blockMeshPtr_(block),
    dict_(smootherDict),
    pointCells_(labelListList(blockMeshPtr_->points().size())),
    cellPoints_(labelListList(blockMeshPtr_->cells().size())),
    pointCellsS_(List<std::set<label> >(blockMeshPtr_->points().size())),
    cellNeighbors_(List<std::set<label> >(blockMeshPtr_->cells().size())),
    sumCellQuality_(scalarList(blockMeshPtr_->points().size())),
    cellQuality_(scalarList(blockMeshPtr_->cells().size())),
    pointTopology_(blockMeshTopology(block, smootherDict.subDict("features"))),
    writeIntermediateMesh_(args.optionFound("writeStep"))
{
    Info<< "\nGETMe smoothing with:"  << nl
        << "  - Max iterations            : "
        << readLabel(dict_.lookup("maxGETMeIter")) << nl
        << "  - Tranformation parameter   : "
        << readScalar(dict_.lookup("elemTransformParameter")) << nl
        << "  - Improvement tolerance     : "
        << readScalar(dict_.lookup("improvementTolerance")) << nl
        << "  - Max ineffective iteration : "
        << readScalar(dict_.lookup("maxIneffectiveIteration")) << nl
        << "  - Mean relaxation table     : "
        << readList<scalar>(dict_.lookup("qMeanRelaxationTable")) << nl
        << "  - Min relaxation table      : "
        << readList<scalar>(dict_.lookup("qMinRelaxationTable")) << nl << nl;

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
}

Foam::blockMeshSmoother::~blockMeshSmoother()
{
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::blockMeshSmoother::meshMeanRatio()
{
    // Reset sumCellQuality_
    sumCellQuality_ = scalarList(blockMeshPtr_->points().size(), 0.0);
    minQuality_ = 1.0;
    meanQuality_ = 0.0;
    forAll (blockMeshPtr_->cells(), cellI)
    {
        pointField H
        (
            blockMeshPtr_->cells()[cellI].points(blockMeshPtr_->points())
        );
        cellQuality_[cellI] = cellSmoother(H).meanRatio();

        if (cellQuality_[cellI] < minQuality_)
        {
            minQuality_ = cellQuality_[cellI]; // Store the new qMin
        }
        meanQuality_ += cellQuality_[cellI];

        forAll(cellPoints_[cellI], pointI)
        {
            sumCellQuality_[cellPoints_[cellI][pointI]] += cellQuality_[cellI];
        }
    }
    meanQuality_ /= blockMeshPtr_->cells().size();
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
            pointField H
            (
                blockMeshPtr_->cells()[cellI].points(blockMeshPtr_->points())
            );

            const pointField Hr(cellSmoother(H).geometricTranform(eTP));

            // Add Transformed Element Nodes And Weights
            forAll (cellPoints_[cellI], pointI)
            {
                // compute the associated weight
                const scalar wji
                (
                    std::sqrt
                    (
                        sumCellQuality_[cellPoints_[cellI][pointI]] /
                        ((cellNeighbors_[cellI].size()*cellQuality_[cellI]) +
                        VSMALL)
                    )
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
    const scalarList &rT,
    labelList &relaxedCells
)
{
    // Relaxed nodes
    pointField pip(blockMeshPtr_->points());

    // Relaxation level
    labelList nR(blockMeshPtr_->points().size(), 0);

    while (!tn.empty())
    {
        // Relax all the points marked for move
        for
        (
            std::set<label>::iterator ptI = tn.begin();
            ptI != tn.end();
            ++ptI
        )
        {
            const scalar r(rT[nR[*ptI]]);

            pip[*ptI] = pointTopology_.optimalPoint
            (
                *ptI,
                (1.0 - r)*blockMeshPtr_->points()[*ptI] + r*pi[*ptI]
            );
        }

        // compute quality with modified mesh
        scalarList cq(blockMeshPtr_->cells().size());
        forAll (blockMeshPtr_->cells(), cellI)
        {
            pointField Ht(cellPoints_[cellI].size());
            forAll (cellPoints_[cellI], ptI)
            {
                Ht[ptI] = pip[cellPoints_[cellI][ptI]];
            }
            cq[cellI] = cellSmoother(Ht).meanRatio();

        }

        // Test new cells
        tn.clear();
        forAll (blockMeshPtr_->cells(), cellI)
        {
            if(cq[cellI] < VSMALL)
            { // Set point to reset
                forAll (cellPoints_[cellI], ptI)
                { // Insert cell point in point to move
                    tn.insert(cellPoints_[cellI][ptI]);
                }
            }
        }

        // Reduce the relaxation factor for each point to correct
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

        // Count number of relaxed cell for display
        if (reversedCells != 0)
        {
            relaxedCells.append(reversedCells);
        }
    }

    return pip;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockMeshSmoother::smoothing(const argList &args)
{
    const scalar eTP(readScalar(dict_.lookup("elemTransformParameter")));

    meshMeanRatio();

    Info<< "\nAnalysis before start: mean quality = " << meanQuality_
        << ",\t min quality = " << minQuality_ << endl;

    label nbIterations(0), noMinImproveCounter(0), cycle(0), prevCycle(cycle);

    scalar targetQual(1.0);
    scalar previousMean(meanQuality_), previousMin(minQuality_);

    scalarList rT = readList<scalar>(dict_.lookup("qMeanRelaxationTable"));
    const scalar improveTol(readScalar(dict_.lookup("improvementTolerance")));

    const label maxNoMinImproveCounter
    (
        readLabel(dict_.lookup("maxIneffectiveIteration"))
    );

    bool converged(false);
    while
    (
        nbIterations < readLabel(dict_.lookup("maxGETMeIter")) &&
        !converged
    )
    {
        if (writeIntermediateMesh_)
        {
            word defaultFacesName = "defaultFaces";
            word defaultFacesType = "patch";

            Time runTime(Foam::Time::controlDictName,  args);
            runTime.setTime(scalar(nbIterations), nbIterations);
            runTime.setDeltaT(scalar(1.0));
            polyMesh mesh
            (
                IOobject
                (
                    "",
                    runTime.caseName(),
                    runTime
                ),
                xferCopy(blockMeshPtr_->points()),    // could we re-use space?
                blockMeshPtr_->cells(),
                blockMeshPtr_->patches(),
                blockMeshPtr_->patchNames(),
                blockMeshPtr_->patchDicts(),
                defaultFacesName,
                defaultFacesType
            );

            // Set the precision of the points data to 10
            IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

            mesh.removeFiles();
            if (!mesh.write())
            {
                FatalErrorIn(args.executable())
                    << "Failed writing polyMesh."
                    << exit(FatalError);
            }

            // TODO add volScalarField quality
        }

        // Storage of Temporary Nodes And Weights sum
        scalarList wj(blockMeshPtr_->points().size(), 0.0);
        std::set<label> tp; // List of transformed points

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

        addUntransformedElementNodesAndWeights(pi, wj, tp, targetQual);

        computeNewNodes(pi, wj, tp);

        labelList nbRelCel;

        // Update the mesh with new points
        blockMeshPtr_->setPoints(iterativeNodeRelaxation(pi, tp, rT, nbRelCel));

        // Compute new min and avg quality
        meshMeanRatio();

        const scalar qMeanImprove(meanQuality_ - previousMean);
        const scalar qMinImprove(minQuality_ - previousMin);

        prevCycle = cycle;

        if (cycle == 0 && qMeanImprove < improveTol && nbIterations != 0)
        { // No mean evolution
            cycle = 2;
            rT = readList<scalar>(dict_.lookup("qMinRelaxationTable"));
        }

        if (cycle == 1)
        { // Min cycle running
            if (qMinImprove < VSMALL)
            { // No improvement
                ++noMinImproveCounter;
            }

            if(noMinImproveCounter > maxNoMinImproveCounter)
            { // More than 5 time with no evolution
                cycle = 2;
            }
        }

        if (cycle == 2)
        { // Min cycle start
            if (qMinImprove < improveTol && prevCycle != 0)
            {
                converged = true;
            }
            else
            {
                cycle = 1;
                noMinImproveCounter = 0;

                scalarList cqs(cellQuality_);
                std::sort(cqs.begin(), cqs.end());

                targetQual = cqs[blockMeshPtr_->cells().size()*0.5];
            }
        }

        switch (prevCycle)
        {
        case 0:
        {
            Info<< "  mean cycle ";
            break;
        }
        case 1:
        {
            Info<< "  min cycle  ";
            break;
        }
        case 2:
        {
            Info<< "  converged  ";
            break;
        }
        }

        Info<< "iter " << nbIterations << "\t mean quality = "
            << meanQuality_ << "\t min quality = " << minQuality_ << endl;

        previousMean = meanQuality_;
        previousMin = minQuality_;

        ++nbIterations;
    }
    if (nbIterations != readLabel(dict_.lookup("maxGETMeIter")))
    {
        Info<< "GETMe smoothing converged in " << nbIterations - 1
            << " iterations" << endl;
    }
    else
    {
        Info<< "GETMe smoothing not converged!" << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
