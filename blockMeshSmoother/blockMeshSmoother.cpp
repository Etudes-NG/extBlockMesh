#include "blockMeshSmoother.h"

#include "Time.H"
#include "cellSmoother.h"

#include <algorithm>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::blockMeshSmoother::blockMeshSmoother
(
    blockMesh &block,
    dictionary &smootherDict,
    const argList &args
)
:
    blockMeshPtr_(&block),
    dict_(smootherDict),
    pointCells_(labelListList(blockMeshPtr_->points().size())),
    cellPoints_(labelListList(blockMeshPtr_->cells().size())),
    pointCellsS_(List<std::set<label> >(blockMeshPtr_->points().size())),
    cellNeighbors_(List<std::set<label> >(blockMeshPtr_->cells().size())),
    sumCellQuality_(scalarList(blockMeshPtr_->points().size())),
    cellQuality_(scalarList(blockMeshPtr_->cells().size())),
    writeIntermediateMesh_(args.optionFound("writeStep")),
    fixBoundary_(args.optionFound("fixBoundary"))
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

    std::set<label> allPoints;
    label val(0);
    forAll (blockMeshPtr_->points(), pointI)
    {
        allPoints.insert(allPoints.end(),val);
        ++val;
    }

    // Create a set of fixed points
    std::set<label> fixedPoints;

    if (fixBoundary_)
    {
        // Search edges
        List<std::set<std::set<label> > > patchFacesPoints
        (
            blockMeshPtr_->patches().size()
        );

        //       points labels              patch  face
        std::map<std::set<label>, std::pair<label, label> > bndFacesLoc;

        //patch  face  points labels
        List<List<std::set<label> > > bndFacesName
        (
            blockMeshPtr_->patches().size()
        );

        forAll (blockMeshPtr_->patches(), patchI)
        {
            bndFacesName[patchI].resize
            (
                blockMeshPtr_->patches()[patchI].size()
            );

            std::set<std::set<label> > facesPointsSet;
            forAll (blockMeshPtr_->patches()[patchI], faceI)
            {
                const labelList facePoints
                (
                    blockMeshPtr_->patches()[patchI][faceI].pointsLabel()
                );
                std::set<label> facePts;
                forAll (facePoints, pointI)
                {
                    facePts.insert(facePoints[pointI]);
                }
                facesPointsSet.insert(facePts);
                bndFacesLoc.insert
                (
                    std::pair<std::set<label>, std::pair<label, label> >
                    (
                        facePts,
                        std::pair<label, label>(patchI, faceI)
                    )
                );
                bndFacesName[patchI][faceI] = facePts;
            }
            patchFacesPoints[patchI] = facesPointsSet;
        }
        //                 patch  face    neiboor set
        List<List<std::set<std::set<label> > > > bndFaceNeiboor
        (
            blockMeshPtr_->patches().size()
        );
        List<List<std::set<std::set<label> > > > bndFaceConnectEdges
        (
            blockMeshPtr_->patches().size()
        );

        forAll (blockMeshPtr_->patches(), patchI)
        {
            bndFaceNeiboor[patchI].resize
            (
                blockMeshPtr_->patches()[patchI].size()
            );
            bndFaceConnectEdges[patchI].resize
            (
                blockMeshPtr_->patches()[patchI].size()
            );

            forAll (blockMeshPtr_->patches()[patchI], faceI)
            { // for all faces of this patch
                forAll (blockMeshPtr_->patches()[patchI], faceJ)
                { // Test all faces of this patch
                    std::set<label> commonPts;
                    std::set_intersection
                    (
                        bndFacesName[patchI][faceI].begin(),
                        bndFacesName[patchI][faceI].end(),
                        bndFacesName[patchI][faceJ].begin(),
                        bndFacesName[patchI][faceJ].end(),
                        std::inserter(commonPts, commonPts.begin())
                    );

                    if (commonPts.size() == 2)
                    { // Faces share an edge
                        bndFaceNeiboor[patchI][faceI].insert(bndFacesName[patchI][faceJ]);
                        bndFaceConnectEdges[patchI][faceI].insert(commonPts);
                    }
                }
            }
        }

        List<std::set<std::set<std::set<label> > > > patchEdges
        (
            blockMeshPtr_->patches().size()
        );

        forAll (blockMeshPtr_->patches(), patchI)
        {
            forAll (blockMeshPtr_->patches()[patchI], faceI)
            {
                if(bndFaceNeiboor[patchI][faceI].size() == 3)
                { // Face as 3 edges attached with this patch
                    std::set<std::set<label> > faceEdge;
                    const labelList ptLabel
                    (
                        blockMeshPtr_->patches()[patchI][faceI].pointsLabel()
                    );
                    forAll (ptLabel, ptI)
                    {
                        std::set<label> edge;
                        edge.insert(ptLabel[ptI]);
                        edge.insert(ptLabel[(ptI + 1) % 4]);
                        faceEdge.insert(edge);
                    }

                    Info<< "Set from pointsLabels   ";
                    for
                    (
                        std::set<std::set<label> >::iterator iter1 = faceEdge.begin();
                        iter1 != faceEdge.end();
                        ++iter1
                    )
                    {
                        for
                        (
                            std::set<label>::iterator iter2 = iter1->begin();
                            iter2 != iter1->end();
                            ++iter2
                        )
                        {
                            Info<< (*iter2) << " ";
                        }
                        Info<< " - ";
                    }
                    Info<< nl;

                    Info<< "Set from bndFaceNeiboor ";
                    for
                    (
                        std::set<std::set<label> >::iterator iter1 = bndFaceConnectEdges[patchI][faceI].begin();
                        iter1 != bndFaceConnectEdges[patchI][faceI].end();
                        ++iter1
                    )
                    {
                        for
                        (
                            std::set<label>::iterator iter2 = iter1->begin();
                            iter2 != iter1->end();
                            ++iter2
                        )
                        {
                            Info<< (*iter2) << " ";
                        }
                        Info<< " - ";
                    }
                    Info<< nl;

                    std::set<std::set<label> > bndEdge;
                    std::set_difference
                    (
                        faceEdge.begin(),
                        faceEdge.end(),
                        bndFaceNeiboor[patchI][faceI].begin(),
                        bndFaceNeiboor[patchI][faceI].end(),
                        std::inserter(bndEdge, bndEdge.begin())
                    );

                    Info<< label(bndEdge.size()) << " - " << label(faceEdge.size()) << endl;
                    patchEdges[patchI].insert(bndEdge);
                }
            }
        }
        forAll (blockMeshPtr_->patches(), patchI)
        {
            Info<< "Patch " << patchI;
            label edg(0);
            for
            (
                std::set<std::set<std::set<label> > >::iterator
                    edgeI = patchEdges[patchI].begin();
                edgeI != patchEdges[patchI].end();
                ++edgeI
            )
            { // for each edge
                Info<< " edge " << edg;
                for
                (
                    std::set<std::set<label> >::iterator ptI = edgeI->begin();
                    ptI != edgeI->end();
                    ++ptI
                )
                {
//                    Info<<;
                }
                ++edg;
            }
            Info<< endl;
        }

        forAll (blockMeshPtr_->patches(), patchI)
        {
            forAll (blockMeshPtr_->patches()[patchI], faceI)
            {
                const labelList facePoints
                (
                    blockMeshPtr_->patches()[patchI][faceI].pointsLabel()
                );
                forAll (facePoints, pointI)
                {
                    fixedPoints.insert(facePoints[pointI]);

                    std::set<label> triangle1;
                    triangle1.insert(facePoints[(pointI + 1) % 4]);
                    triangle1.insert(facePoints[(pointI + 2) % 4]);

                    std::set<label> triangle2;
                    triangle2.insert(facePoints[(pointI + 2) % 4]);
                    triangle2.insert(facePoints[(pointI + 3) % 4]);

                    std::map<label, std::set<std::set<label> > >::iterator
                         iter = bndPtTri_.find(facePoints[pointI]);

                    if (iter != bndPtTri_.end())
                    {
                        iter->second.insert(triangle1);
                        iter->second.insert(triangle2);
                    }
                    else
                    {
                        std::set<std::set<label> > triangleList;
                        triangleList.insert(triangle1);
                        triangleList.insert(triangle2);
                        bndPtTri_[facePoints[pointI]] = triangleList;
                    }
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

        // fixed point boundary faces topology
        for
        (
            std::set<label>::iterator fixPtI = fixedPoints.begin();
            fixPtI != fixedPoints.end();
            ++fixPtI
        )
        {
//            blockMeshPtr_->points()[*fixPtI].
        }
    }
    else
    {
        mobilPoints_ = allPoints;
    }
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
        cellQuality_[cellI] =
        cellSmoother
        (
            blockMeshPtr_->cells()[cellI].points
            (
             blockMeshPtr_->points()
            )
        ).meanRatio();

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
            pip[*ptI] = (1.0 - r)*blockMeshPtr_->points()[*ptI] + r*pi[*ptI];
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

        // Remove fixed points from point to move if any
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
        readLabel
        (
            dict_.lookup("maxIneffectiveIteration")
        )
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

        labelList nbRelCel;

        // Update the mesh with new points
        blockMeshPtr_->setPoints(iterativeNodeRelaxation(pi, tn, rT, nbRelCel));

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
