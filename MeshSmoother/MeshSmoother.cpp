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

#include "MeshSmoother.h"

#include "triSurface.H"
#include "polyPatch.H"
#include "Time.H"
#include "IOmanip.H"
#include "blockMesh.H"

#include "SmootherPoint.h"
#include "SmootherCell.h"
#include "SmootherControl.h"
#include "SmootherParameter.h"
#include "SmootherBoundary.h"

#include <algorithm>
#include <cmath>

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void MeshSmoother::analyseMeshQuality()
{
    forAll(_cell, cellI)
    {
        _cell[cellI]->computeQuality();
    }
}

void Foam::MeshSmoother::analyseMeshQuality(const labelHashSet &cell)
{
    forAllConstIter(labelHashSet , cell, cellI)
    {
        _cell[cellI.key()]->computeQuality();
    }
}

void Foam::MeshSmoother::qualityStats()
{
    scalar minQuality = 1.0;
    scalar meanQuality = 0.0;

    scalarList pQSum(_polyMesh->nPoints(), 0.0);
    forAll(_polyMesh->cells(), cellI)
    {
        const scalar& cQ = _cell[cellI]->quality();
        if (cQ < minQuality)
        {
            minQuality = cQ;
        }

        meanQuality += cQ;

        forAll(_polyMesh->cellPoints()[cellI], pointI)
        {
            pQSum[_polyMesh->cellPoints()[cellI][pointI]] += cQ;
        }
    }
    forAll(_polyMesh->points(), ptI)
    {
        _bnd->pt(ptI)->setQuality
        (
            pQSum[ptI]/_polyMesh->pointCells()[ptI].size()
        );
    }

    meanQuality /= _polyMesh->nCells();

    _param->setMinQual(minQuality);
    _param->setMeanQual(meanQuality);
}

Foam::labelHashSet Foam::MeshSmoother::addTransformedElementNodeWeight()
{
    labelHashSet transformedPoints;
    forAll(_polyMesh->cells(), cellI)
    {
        if (_cell[cellI]->quality() <= _param->transformationTreshold())
        {
            const pointField newCellPoints = _cell[cellI]->geometricTransform();
            const cellShape& cS = _polyMesh->cellShapes()[cellI];
            const scalar& cQ = _cell[cellI]->quality();

            if (cQ < VSMALL)
            {
                FatalErrorIn("Foam::MeshSmoother::weightingFactor()")
                    << "Quality of cell " << cellI << " is null" << nl
                    << exit(FatalError);
            }

            forAll(cS, pointI)
            {
                SmootherPoint* pt = _bnd->pt(cS[pointI]);

                // compute the associated weight
                const label nNei = _polyMesh->pointCells(cS[pointI]).size();
                const scalar weight = std::sqrt(pt->avgQual()/(nNei*cQ));

                // add the weight to temporary weighted sum
                pt->addWeight(weight, newCellPoints[pointI]);

                // add point to set of tranformed points
                transformedPoints.insert(cS[pointI]);
            }
        }
    }
    return transformedPoints;
}

void Foam::MeshSmoother::addUnTransformedElementNodeWeight(labelHashSet &tp)
{
    // Add Untransformed Element Nodes And Weights
    forAll (_polyMesh->cells(), cellI)
    {
        if (untransformedAndhavePointTransformed(cellI, tp))
        {
            const cellShape& cS = _polyMesh->cellShapes()[cellI];
            const scalar& cQ = _cell[cellI]->quality();

            if (cQ < VSMALL)
            {
                FatalErrorIn("Foam::MeshSmoother::weightingFactor()")
                    << "Quality of cell " << cellI << " is null" << nl
                    << exit(FatalError);
            }

            forAll (cS, pointI)
            {
                SmootherPoint* pt = _bnd->pt(cS[pointI]);

                // compute the associated weight
                const label nNei = _polyMesh->pointCells(cS[pointI]).size();
                const scalar weight = std::sqrt(pt->avgQual()/(nNei*cQ));

                // add the weight to temporary weighted sum
                pt->addWeight(weight);
            }
        }
    }
}

bool Foam::MeshSmoother::untransformedAndhavePointTransformed
(
    const label cellI,
    const labelHashSet& tp
)
{
    if (_cell[cellI]->quality() > _param->transformationTreshold())
    {
        const cellShape& cS = _polyMesh->cellShapes()[cellI];
        forAll(cS, ptI)
        {
            if (tp.find(cS[ptI]) != tp.end())
            { // Have point transformed

                return true;
            }
        }
    }

    return false;
}

void Foam::MeshSmoother::iterativeNodeRelaxation
(
    labelHashSet &tP,
    const scalarList &r
)
{
    // Reset relaxation level
    forAll(_polyMesh->points(), ptI)
    {
        _bnd->pt(ptI)->resetRelaxationLevel();
    }
    _param->setNbMovedPoints(tP.size());

    label nbRelax = 0;
    while (!tP.empty())
    {
        ++nbRelax;
        // Relax all the points marked for move
        labelHashSet modifiedCells;
        forAllConstIter(labelHashSet, tP, ptI)
        {
            _bnd->pt(ptI.key())->relaxPoint(r);

            forAll(_polyMesh->pointCells()[ptI.key()], cellI)
            {
                modifiedCells.insert(_polyMesh->pointCells()[ptI.key()][cellI]);
            }
        }

        // compute quality with relaxed points
        tP.clearStorage();
        analyseMeshQuality(modifiedCells);
        forAllConstIter(labelHashSet , modifiedCells, cellI)
        {
            if(_cell[cellI.key()]->quality() < VSMALL)
            {
                tP.insert(_polyMesh->cellPoints()[cellI.key()]);
            }
        }

        // Increase the relaxation level for invalid points
        forAllIter(labelHashSet, tP, ptI)
        {
            _bnd->pt(ptI.key())->addRelaxLevel(r);
        }
    }
    _param->setNbRelaxations(nbRelax);
}

bool Foam::MeshSmoother::runIteration()
{
    _param->resetUpdateTime();

    // Store mean and min quality before iteration
    const scalar minQ = _param->minQual();
    const scalar meanQ = _param->meanQual();

    if (!_bnd->unSnapedPoints().empty())
    {
        snapSmoothing();
    }
    else
    {
        GETMeSmoothing();
    }

    // Compute new min and avg quality
    qualityStats();
    _param->printStatus(_bnd->unSnapedPoints().size());

    const bool asUnSnaped = _bnd->unSnapedPoints().empty();
    return _param->setSmoothCycle(meanQ, minQ, asUnSnaped, this);
}

pointField MeshSmoother::getMovedPoints() const
{
    pointField pt(_polyMesh->nPoints());

    forAll(pt, ptI)
    {
        pt[ptI] = _bnd->pt(ptI)->getRelaxedPoint();
    }
    return pt;
}

void Foam::MeshSmoother::writeMesh
(
    const fvMesh& meshFv,
    volScalarField& meshQuality
) const
{
    forAll(meshQuality, cellI)
    {
        meshQuality[cellI] = _cell[cellI]->quality();
    }

    if (!meshFv.write())
    {
        FatalErrorIn("Foam::MeshSmoother::writeMesh()")
            << "Failed writing fvMesh."
            << exit(FatalError);
    }

    if (!_polyMesh->write())
    {
        FatalErrorIn("Foam::MeshSmoother::writeMesh()")
            << "Failed writing polyMesh."
            << exit(FatalError);
    }
}

void MeshSmoother::GETMeSmoothing()
{
    //-------------------------------------------------------------------------

    // Reset all points
    forAll(_polyMesh->points(), ptI)
    {
        _bnd->pt(ptI)->laplaceReset();
    }

    // LaplaceSmooth boundary points
    labelHashSet snapPoints = _bnd->featuresPoints();
    forAllConstIter(labelHashSet, snapPoints, ptI)
    {
        _bnd->pt(ptI.key())->featLaplaceSmooth();
    }
    iterativeNodeRelaxation(snapPoints, _ctrl->snapRelaxTable());

    // Snap boundary points
    snapPoints = _bnd->featuresPoints();
    forAllConstIter(labelHashSet, snapPoints, ptI)
    {
        _bnd->pt(ptI.key())->snap();
    }
    iterativeNodeRelaxation(snapPoints, _ctrl->snapRelaxTable());

    //-------------------------------------------------------------------------

    forAll(_polyMesh->points(), ptI)
    {
        _bnd->pt(ptI)->GETMeReset();
    }

    labelHashSet transformedPoints = addTransformedElementNodeWeight();

    if (transformedPoints.size() != _polyMesh->nPoints())
    {
        addUnTransformedElementNodeWeight(transformedPoints);
    }

    // Compute new point
    forAllConstIter(labelHashSet, transformedPoints, ptI)
    {
        _bnd->pt(ptI.key())->GETMeSmooth();
    }

    iterativeNodeRelaxation(transformedPoints, _param->relaxationTable());

//    _bnd->writeAllSurfaces(_param->getIterNb());
}

void MeshSmoother::snapSmoothing()
{
    // Reset all points
    forAll(_polyMesh->points(), ptI)
    {
        _bnd->pt(ptI)->laplaceReset();
    }

    // LaplaceSmooth interior points
    labelHashSet laplacePoints = _bnd->interiorPoints();
    forAllConstIter(labelHashSet, laplacePoints, ptI)
    {
        _bnd->pt(ptI.key())->laplaceSmooth();
    }
    iterativeNodeRelaxation(laplacePoints, _ctrl->snapRelaxTable());

    // Snap boundary points
    labelHashSet snapPoints = _bnd->featuresPoints();
    forAllConstIter(labelHashSet, snapPoints, ptI)
    {
        _bnd->pt(ptI.key())->snap();
    }
    iterativeNodeRelaxation(snapPoints, _ctrl->snapRelaxTable());

    // Remove points from unsnaped point list if snaped
    forAll(_polyMesh->points(), ptI)
    {
        _bnd->pt(ptI)->needSnap();
    }

    // LaplaceSmooth boundary points
    snapPoints = _bnd->featuresPoints();
    forAllConstIter(labelHashSet, snapPoints, ptI)
    {
        _bnd->pt(ptI.key())->featLaplaceSmooth();
    }
    iterativeNodeRelaxation(snapPoints, _ctrl->snapRelaxTable());
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MeshSmoother::MeshSmoother
(
    polyMesh *mesh,
    dictionary *smootherDict,
    blockMesh *blocks
)
:
    _polyMesh(mesh),
    _blocks(blocks)
{
    scalar time = _polyMesh->time().elapsedCpuTime();

    _ctrl = new SmootherControl(smootherDict);
    _param = new SmootherParameter(_ctrl, _polyMesh);
    dictionary& snapDict = smootherDict->subDict("snapControls");
    _bnd = new SmootherBoundary(snapDict, _polyMesh);
    _cell = List<SmootherCell*>(_polyMesh->nCells());

    SmootherPoint dummyPoint;
    dummyPoint.setStaticItems(_bnd, _param, _polyMesh);

    forAll(_cell, cellI)
    {
        _cell[cellI] = new SmootherCell(_polyMesh->cellShapes()[cellI]);
    }
    _cell[0]->setStaticItems(_bnd,_ctrl->transformationParameter());

    // Analyse initial quality
    analyseMeshQuality();
    qualityStats();

    //snapFeatures();

    Info<< "Smoother initialized in "
        << _polyMesh->time().elapsedCpuTime() - time << " s" << nl << nl
        << "Smooth the mesh" << nl << nl;

    _param->resetUpdateTime();
    _param->printHeaders();
    _param->setSmoothCycle(!_bnd->unSnapedPoints().empty());
    _param->printStatus(_bnd->unSnapedPoints().size());
    _param->setIterNb();
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::MeshSmoother::~MeshSmoother()
{
    forAll(_cell, cellI)
    {
        delete _cell[cellI];
    }

    delete _param;
    delete _bnd;
    delete _ctrl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MeshSmoother::update()
{
    while(runIteration()){}

    // Update the mesh with new points
    _polyMesh->movePoints(getMovedPoints());

    _param->printStats();
}

void Foam::MeshSmoother::updateAndWrite
(
    word& regionName,
    word& defaultFacesName,
    word& defaultFacesType,
    Time& runTime
)
{
    fvMesh meshFv
    (
        IOobject(regionName, runTime.constant(), runTime),
        xferCopy(_polyMesh->points()),
        _blocks->cells(),
        _blocks->patches(),
        _blocks->patchNames(),
        _blocks->patchDicts(),
        defaultFacesName,
        defaultFacesType
    );

    volScalarField meshQuality
    (
        IOobject
        (
            "meshQuality",
            runTime.timeName(),
            meshFv,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        meshFv,
        0.0
    );

    writeMesh(meshFv, meshQuality);

    while(runIteration())
    {
        // Update the mesh with new points
        _polyMesh->movePoints(getMovedPoints());

        ++runTime;
        writeMesh(meshFv, meshQuality);
    }

    _param->printStats();
}

Foam::scalar Foam::MeshSmoother::getTransformationTreshold() const
{
    scalarList cqs(_polyMesh->nCells());
    forAll(_polyMesh->cells(), cellI)
    {
        cqs[cellI] = _cell[cellI]->quality();
    }
    std::sort(cqs.begin(), cqs.end());

    return cqs[std::floor(_polyMesh->nCells()*_ctrl->ratioForMin())];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
