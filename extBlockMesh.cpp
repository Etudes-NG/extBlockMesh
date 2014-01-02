/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    blockMesh

Description
    A multi-block mesh generator.

    Uses the block mesh description found in
    \a constant/polyMesh/blockMeshDict
    (or \a constant/\<region\>/polyMesh/blockMeshDict).

Usage

    - blockMesh [OPTION]

    \param -blockTopology \n
    Write the topology as a set of edges in OBJ format.

    \param -region \<name\> \n
    Specify an alternative mesh region.

    \param -dict \<filename\> \n
    Specify alternative dictionary for the block mesh description.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "IOdictionary.H"
#include "IOPtrList.H"

#include "attachPolyTopoChanger.H"
#include "emptyPolyPatch.H"
#include "cellSet.H"

#include "argList.H"
#include "OSspecific.H"
#include "OFstream.H"

#include "Pair.H"
#include "slidingInterface.H"

#include "blockMeshSmoother.h"
//-----------------------------------------


#include <set>
#include <vector>
#include <algorithm>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar meanRatio(const pointField &H)
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
    forAll (H, ptI)
    {
        const point p1(H[v1[ptI]] - H[ptI]);
        const point p2(H[v2[ptI]] - H[ptI]);
        const point p3(H[v3[ptI]] - H[ptI]);
        const Tensor<scalar> mA(p1, p2, p3);

        const scalar sigma(det(mA));

        if (sigma > VSMALL)
        {
            cQa += 3*std::pow(sigma, 2.0/3.0)/magSqr(mA);
        }
    }
    return cQa/8.0;
}

pointField geometricTransformation
(
    const pointField &H,
    const scalar &cor
)
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

    // Labels for edge lengh
    labelList vd1(12),  vd2(12);
    vd1[0] = 0;	vd1[1] = 1;	vd1[2] = 2;	vd1[3] = 3;
    vd1[4] = 0;	vd1[5] = 1;	vd1[6] = 2;	vd1[7] = 3;
    vd1[8] = 4;	vd1[9] = 5;	vd1[10]= 6; vd1[11]= 7;

    vd2[0] = 1;	vd2[1] = 2;	vd2[2] = 3;	vd2[3] = 0;
    vd2[4] = 4;	vd2[5] = 5;	vd2[6] = 6;	vd2[7] = 7;
    vd2[8] = 5;	vd2[9] = 6;	vd2[10]= 7;	vd2[11]= 4;

    // Compute dual octahedron
    pointField oct(6);
    forAll (vb1, octPtI)
    {
        oct[octPtI] =
        (
            H[vb1[octPtI]] + H[vb2[octPtI]] +
            H[vb3[octPtI]] + H[vb4[octPtI]]
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
    forAll (H, ptI)
    {
        Hp[ptI] = octC[ptI] + cor/std::sqrt(mag(octN[ptI]))*octN[ptI];
    }

    // Scaling of new points
    // Centroid of cell
    point c(0, 0, 0);
    forAll (H, ptI)
    {
        c += Hp[ptI];
    }
    c /= 8;

    // Scaling factor (keeping avg edge lenght)
    point mh1(0, 0, 0), mh2(0, 0, 0);
    forAll (vd1, edgeI)
    {
        mh1 += H[vd1[edgeI]] - H[vd2[edgeI]];
        mh2 += Hp[vd1[edgeI]] - Hp[vd2[edgeI]];
    }
    mh1 /= 12;
    mh2 /= 12;
    const scalar scalingfact(mag(mh1)/mag(mh2));

    pointField C(8, c);

    // Return relaxation of new points
    return C + scalingfact*(Hp - C);
}

pointField pointsData(const labelList &labels, const pointField &points)
{
    pointField out(labels.size());
    forAll (labels, pointI)
    {
        out[pointI] = points[labels[pointI]];
    }
    return out;
}

std::set<label> pointsToRevert
(
    const pointField &pts,
    const labelListList &cp,
    const scalarList &cq
)
{
    std::set<label> pointsToRevert;
    forAll (cq, cellI)
    {
        if (meanRatio(pointsData(cp[cellI], pts)) < VSMALL /*cq[cellI]*/)
        {
            forAll (cp[cellI], pointI)
            {
                pointsToRevert.insert(cp[cellI][pointI]);
            }
        }
    }
    return pointsToRevert;
}

scalarList meshMeanRatio
(
    scalar &min,
    scalar &avg,
    scalarList &pwi,
    const blockMesh &blocks,
    const labelListList &cp
)
{
    scalarList cq(blocks.cells().size());
    min = 1.0;
    avg = 0;
    forAll (blocks.cells(), cellI)
    {
        cq[cellI] = meanRatio(blocks.cells()[cellI].points(blocks.points()));

        if (cq[cellI] < min)
        {
            min = cq[cellI]; // Store the new qMin
        }
        avg += cq[cellI];

        forAll(cp[cellI], pointI)
        {
            pwi[cp[cellI][pointI]] += cq[cellI];
        }
    }
    avg /= blocks.cells().size();
    return cq;
}

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption
    (
        "blockTopology",
        "write block edges and centres as .obj files"
    );
    argList::addOption
    (
        "dict",
        "file",
        "specify alternative dictionary for the blockMesh description"
    );

#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"

    const word dictName("blockMeshDict");

    word regionName;
    fileName polyMeshDir;

    if (args.optionReadIfPresent("region", regionName, polyMesh::defaultRegion))
    {
        // constant/<region>/polyMesh/blockMeshDict
        polyMeshDir = regionName/polyMesh::meshSubDir;

        Info<< nl << "Generating mesh for region " << regionName << endl;
    }
    else
    {
        // constant/polyMesh/blockMeshDict
        polyMeshDir = polyMesh::meshSubDir;
    }

    IOobject meshDictIO
    (
        dictName,
        runTime.constant(),
        polyMeshDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (args.optionFound("dict"))
    {
        const fileName dictPath = args["dict"];

        meshDictIO = IOobject
        (
            (
                isDir(dictPath)
              ? dictPath/dictName
              : dictPath
            ),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );
    }

    if (!meshDictIO.headerOk())
    {
        FatalErrorIn(args.executable())
            << "Cannot open mesh description file\n    "
            << meshDictIO.objectPath()
            << nl
            << exit(FatalError);
    }

    Info<< "Creating block mesh from\n    "
        << meshDictIO.objectPath() << endl;

    blockMesh::verbose(true);

    IOdictionary meshDict(meshDictIO);
    blockMesh blocks(meshDict, regionName);

    // GETMe adaptive smoothing
    if (meshDict.found("smoother"))
    {
        dictionary smoothDict(meshDict.subDict("smoother"));
        blockMeshSmoother smoother(blocks, smoothDict);
        smoother.smoothing();
    }

    if (args.optionFound("blockTopology"))
    {
        // Write mesh as edges.
        {
            fileName objMeshFile("blockTopology.obj");

            OFstream str(runTime.path()/objMeshFile);

            Info<< nl << "Dumping block structure as Lightwave obj format"
                << " to " << objMeshFile << endl;

            blocks.writeTopology(str);
        }

        // Write centres of blocks
        {
            fileName objCcFile("blockCentres.obj");

            OFstream str(runTime.path()/objCcFile);

            Info<< nl << "Dumping block centres as Lightwave obj format"
                << " to " << objCcFile << endl;

            const polyMesh& topo = blocks.topology();

            const pointField& cellCentres = topo.cellCentres();

            forAll(cellCentres, cellI)
            {
                //point cc = b.blockShape().centre(b.points());
                const point& cc = cellCentres[cellI];

                str << "v " << cc.x() << ' ' << cc.y() << ' ' << cc.z() << nl;
            }
        }

        Info<< nl << "end" << endl;

        return 0;
    }


//    Info<< nl << "Creating polyMesh from blockMesh" << endl;

    word defaultFacesName = "defaultFaces";
    word defaultFacesType = emptyPolyPatch::typeName;
    polyMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.constant(),
            runTime
        ),
        xferCopy(blocks.points()),           // could we re-use space?
        blocks.cells(),
        blocks.patches(),
        blocks.patchNames(),
        blocks.patchDicts(),
        defaultFacesName,
        defaultFacesType
    );

    // Read in a list of dictionaries for the merge patch pairs
    if (meshDict.found("mergePatchPairs"))
    {
        List<Pair<word> > mergePatchPairs
        (
            meshDict.lookup("mergePatchPairs")
        );

#       include "mergePatchPairs.H"
    }
    else
    {
//        Info<< nl << "There are no merge patch pairs edges" << endl;
    }


    // Set any cellZones (note: cell labelling unaffected by above
    // mergePatchPairs)

    label nZones = blocks.numZonedBlocks();

    if (nZones > 0)
    {
        Info<< nl << "Adding cell zones" << endl;

        // Map from zoneName to cellZone index
        HashTable<label> zoneMap(nZones);

        // Cells per zone.
        List<DynamicList<label> > zoneCells(nZones);

        // Running cell counter
        label cellI = 0;

        // Largest zone so far
        label freeZoneI = 0;

        forAll(blocks, blockI)
        {
            const block& b = blocks[blockI];
            const labelListList& blockCells = b.cells();
            const word& zoneName = b.blockDef().zoneName();

            if (zoneName.size())
            {
                HashTable<label>::const_iterator iter = zoneMap.find(zoneName);

                label zoneI;

                if (iter == zoneMap.end())
                {
                    zoneI = freeZoneI++;

                    Info<< "    " << zoneI << '\t' << zoneName << endl;

                    zoneMap.insert(zoneName, zoneI);
                }
                else
                {
                    zoneI = iter();
                }

                forAll(blockCells, i)
                {
                    zoneCells[zoneI].append(cellI++);
                }
            }
            else
            {
                cellI += b.cells().size();
            }
        }


        List<cellZone*> cz(zoneMap.size());

//        Info<< nl << "Writing cell zones as cellSets" << endl;

        forAllConstIter(HashTable<label>, zoneMap, iter)
        {
            label zoneI = iter();

            cz[zoneI] = new cellZone
            (
                iter.key(),
                zoneCells[zoneI].shrink(),
                zoneI,
                mesh.cellZones()
            );

            // Write as cellSet for ease of processing
            cellSet cset(mesh, iter.key(), zoneCells[zoneI].shrink());
            cset.write();
        }

        mesh.pointZones().setSize(0);
        mesh.faceZones().setSize(0);
        mesh.cellZones().setSize(0);
        mesh.addZones(List<pointZone*>(0), List<faceZone*>(0), cz);
    }

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

//    Info<< nl << "Writing polyMesh" << endl;
    mesh.removeFiles();
    if (!mesh.write())
    {
        FatalErrorIn(args.executable())
            << "Failed writing polyMesh."
            << exit(FatalError);
    }


    //
    // write some information
    //
    {
        const polyPatchList& patches = mesh.boundaryMesh();

//        Info<< "----------------" << nl
//            << "Mesh Information" << nl
//            << "----------------" << nl
//            << "  " << "boundingBox: " << boundBox(mesh.points()) << nl
//            << "  " << "nPoints: " << mesh.nPoints() << nl
//            << "  " << "nCells: " << mesh.nCells() << nl
//            << "  " << "nFaces: " << mesh.nFaces() << nl
//            << "  " << "nInternalFaces: " << mesh.nInternalFaces() << nl;

//        Info<< "----------------" << nl
//            << "Patches" << nl
//            << "----------------" << nl;

        forAll(patches, patchI)
        {
            const polyPatch& p = patches[patchI];

//            Info<< "  " << "patch " << patchI
//                << " (start: " << p.start()
//                << " size: " << p.size()
//                << ") name: " << p.name()
//                << nl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
