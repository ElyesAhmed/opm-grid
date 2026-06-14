/*
  Copyright 2026 SINTEF Digital.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/grid/cpgrid/refinement/LeafGridAssembler.hpp>

#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/EntityRep.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <cstdio>
#include <limits>

namespace
{

using Dune::cpgrid::CpGridData;
using Dune::cpgrid::EntityRep;
using Opm::Refinement::BlockRefinement;
using Opm::Refinement::GridStateWriter;

struct SourceRef
{
    int grid;  // 0 = level zero, b+1 = level grid of box b
    int index; // cell or face index in that grid
};

int axisOf(face_tag tag)
{
    switch (tag) {
    case I_FACE: return 0;
    case J_FACE: return 1;
    case K_FACE: return 2;
    default:
        throw std::logic_error("Unexpected face tag in refined level grid.");
    }
}

// In a distributed grid, face_to_cell uses this sentinel for a neighbor cell
// that lives on another rank (not present locally).
constexpr int kRemoteCell = std::numeric_limits<int>::max();

} // anonymous namespace

namespace Opm
{
namespace Refinement
{

std::shared_ptr<CpGridData>
assembleLeafGrid(std::vector<std::shared_ptr<CpGridData>>& storage,
                 const std::vector<BlockRefinement>& requests,
                 Dune::MPIHelper::MPICommunicator comm)
{
    CpGridData& level0 = *storage[0];
    const int numBoxes = static_cast<int>(requests.size());
    const auto& dims0 = level0.logicalCartesianSize();
    const int numCells0 = level0.size(0);
    const int numCorners0 = level0.size(3);

    // ----- Basic level-zero lookups ------------------------------------
    std::vector<int> compressed0(static_cast<std::size_t>(dims0[0])*dims0[1]*dims0[2], -1);
    for (int c = 0; c < numCells0; ++c) {
        compressed0[level0.globalCell()[c]] = c;
    }
    const auto boxOfIJK = [&](int i, int j, int k) {
        for (int b = 0; b < numBoxes; ++b) {
            const auto& req = requests[b];
            if (i >= req.startIJK[0] && i < req.endIJK[0]
                && j >= req.startIJK[1] && j < req.endIJK[1]
                && k >= req.startIJK[2] && k < req.endIJK[2]) {
                return b;
            }
        }
        return -1;
    };
    // Box membership of each level-zero cell (-1 if unrefined).
    std::vector<int> boxOfCell(numCells0, -1);
    for (int c = 0; c < numCells0; ++c) {
        const int cart = level0.globalCell()[c];
        const int i = cart % dims0[0];
        const int j = (cart / dims0[0]) % dims0[1];
        const int k = cart / (dims0[0]*dims0[1]);
        boxOfCell[c] = boxOfIJK(i, j, k);
    }

    auto& cellToPoint0 = GridStateWriter::cellToPoint(level0);
    auto& cellToFace0 = GridStateWriter::cellToFace(level0);
    auto& faceToCell0 = GridStateWriter::faceToCell(level0);
    auto& faceToPoint0 = GridStateWriter::faceToPoint(level0);
    auto& faceTag0 = GridStateWriter::faceTag(level0);
    auto& faceNormals0 = GridStateWriter::faceNormals(level0);
    const int numFaces0 = faceToCell0.size();

    // ----- Per-box data -------------------------------------------------
    struct BoxData
    {
        CpGridData* level;
        std::array<int,3> factors;
        std::array<int,3> refinedDims;
        // children of each level-zero parent, ordered by idxInParent
        std::map<int, std::vector<int>> childrenOfParent;
        // level corner -> level-zero corner (or -1)
        std::vector<int> cornerEquiv;
        // level corner -> leaf pool index
        std::vector<int> cornerToLeaf;
        // level cell/face -> leaf index
        std::vector<int> cellToLeaf;
        std::vector<int> faceToLeaf;
    };
    std::vector<BoxData> boxes(numBoxes);

    for (int b = 0; b < numBoxes; ++b) {
        BoxData& box = boxes[b];
        box.level = storage[b + 1].get();
        box.factors = requests[b].cellsPerDim;
        box.refinedDims = box.level->logicalCartesianSize();

        const auto& childToParent = GridStateWriter::childToParent(*box.level);
        const auto& idxInParent = GridStateWriter::idxInParent(*box.level);
        std::map<int, std::vector<std::pair<int,int>>> tmp; // parent -> (idxInParent, child)
        for (int cell = 0; cell < box.level->size(0); ++cell) {
            tmp[childToParent[cell][1]].emplace_back(idxInParent[cell], cell);
        }
        for (auto& [parent, kids] : tmp) {
            std::sort(kids.begin(), kids.end());
            auto& target = box.childrenOfParent[parent];
            target.reserve(kids.size());
            for (const auto& [pos, cell] : kids) {
                target.push_back(cell);
            }
        }

        // Corner identification with level zero, through the parent cells
        // (per-cell identification stays consistent across faults inside
        // the block, where lattice positions are not unique).
        const auto& [rx, ry, rz] = box.factors;
        box.cornerEquiv.assign(box.level->size(3), -1);
        auto& cellToPointL = GridStateWriter::cellToPoint(*box.level);
        for (int cell = 0; cell < box.level->size(0); ++cell) {
            const int parent = childToParent[cell][1];
            const int pos = idxInParent[cell];
            const int ii = pos % rx;
            const int jj = (pos / rx) % ry;
            const int kk = pos / (rx*ry);
            for (int corner = 0; corner < 8; ++corner) {
                const int di = corner & 1;
                const int dj = (corner >> 1) & 1;
                const int dk = (corner >> 2) & 1;
                if ((ii + di) % rx == 0 && (jj + dj) % ry == 0 && (kk + dk) % rz == 0) {
                    const int pdi = (ii + di) / rx;
                    const int pdj = (jj + dj) / ry;
                    const int pdk = (kk + dk) / rz;
                    const int parentCorner = cellToPoint0[parent][pdi + 2*pdj + 4*pdk];
                    box.cornerEquiv[cellToPointL[cell][corner]] = parentCorner;
                }
            }
        }
    }

    // ----- Corner pool (D1) ----------------------------------------------
    // Pool = all level-zero corners (indices preserved) + new refined corners.
    auto leaf = std::make_shared<CpGridData>(comm, storage);
    auto& leafGeometry = GridStateWriter::geometry(*leaf);
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& leafCorners =
        *(leafGeometry.geomVector(std::integral_constant<int,3>()));
    const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& corners0 =
        *(GridStateWriter::geometry(level0).geomVector(std::integral_constant<int,3>()));
    leafCorners.assign(corners0.begin(), corners0.end());

    std::vector<std::array<int,2>> leafCornerHistory(numCorners0);
    for (int corner = 0; corner < numCorners0; ++corner) {
        leafCornerHistory[corner] = {0, corner};
    }

    // Refined corners shared between touching boxes (edge/corner-sharing)
    // coincide bitwise: both boxes resample the same parent description with
    // equal subdivisions, so the arithmetic is identical. Exact-coordinate
    // matching is therefore correct (not a fragile floating-point heuristic)
    // and merges those corners into one pool entry. Disjoint boxes only
    // coincide on shared boundaries, so a global map over refined corners is
    // safe.
    std::map<std::array<double,3>, int> refinedCornerPool;
    const auto coordKey = [](const Dune::cpgrid::Geometry<0,3>& corner) {
        const auto& c = corner.center();
        return std::array<double,3>{ c[0], c[1], c[2] };
    };

    for (int b = 0; b < numBoxes; ++b) {
        BoxData& box = boxes[b];
        const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& cornersL =
            *(GridStateWriter::geometry(*box.level).geomVector(std::integral_constant<int,3>()));
        box.cornerToLeaf.assign(box.level->size(3), -1);
        std::vector<std::array<int,2>> levelCornerHistory(box.level->size(3), {-1, -1});
        for (int corner = 0; corner < box.level->size(3); ++corner) {
            if (box.cornerEquiv[corner] >= 0) {
                box.cornerToLeaf[corner] = box.cornerEquiv[corner];
                levelCornerHistory[corner] = {0, box.cornerEquiv[corner]};
            }
            else {
                const auto key = coordKey(cornersL[corner]);
                const auto existing = refinedCornerPool.find(key);
                if (existing != refinedCornerPool.end()) {
                    // Shared with an already-processed box: reuse, and keep
                    // the corner's birth identity from that earlier box.
                    box.cornerToLeaf[corner] = existing->second;
                    levelCornerHistory[corner] = leafCornerHistory[existing->second];
                }
                else {
                    const int leafCorner = static_cast<int>(leafCorners.size());
                    box.cornerToLeaf[corner] = leafCorner;
                    leafCorners.push_back(cornersL[corner]);
                    leafCornerHistory.push_back({b + 1, corner});
                    refinedCornerPool.emplace(key, leafCorner);
                }
            }
        }
        GridStateWriter::setCornerHistory(*box.level, std::move(levelCornerHistory));
    }
    const int numLeafCorners = static_cast<int>(leafCorners.size());

    // ----- Leaf cells -----------------------------------------------------
    // Level-zero order, parents replaced by their children (idxInParent order).
    std::vector<SourceRef> leafCells;
    std::vector<int> leafIdxOfCell0(numCells0, -1);
    std::vector<std::vector<int>> leafIdxOfLevelCell(numBoxes);
    for (int b = 0; b < numBoxes; ++b) {
        leafIdxOfLevelCell[b].assign(boxes[b].level->size(0), -1);
    }
    std::vector<std::tuple<int, std::vector<int>>> parentToChildren(numCells0,
                                                                    std::make_tuple(-1, std::vector<int>{}));
    for (int c = 0; c < numCells0; ++c) {
        const int b = boxOfCell[c];
        if (b < 0) {
            leafIdxOfCell0[c] = static_cast<int>(leafCells.size());
            leafCells.push_back({0, c});
        }
        else {
            auto childrenIt = boxes[b].childrenOfParent.find(c);
            if (childrenIt == boxes[b].childrenOfParent.end()) {
                continue; // active parent without active children cannot happen; inactive handled by absence
            }
            for (const int child : childrenIt->second) {
                leafIdxOfLevelCell[b][child] = static_cast<int>(leafCells.size());
                leafCells.push_back({b + 1, child});
            }
            parentToChildren[c] = std::make_tuple(b + 1, childrenIt->second);
        }
    }
    const int numLeafCells = static_cast<int>(leafCells.size());

    // ----- Identify replaced level-zero faces and mosaic pairings ---------
    // A level-zero face is dropped if any adjacent cell is a refined parent.
    // For parent faces toward an unrefined neighbor, the refined boundary
    // faces of the box take over, paired with that neighbor.
    std::vector<bool> faceDropped(numFaces0, false);
    for (int face = 0; face < numFaces0; ++face) {
        const auto row = faceToCell0[EntityRep<1>(face, true)];
        for (int e = 0; e < row.size(); ++e) {
            const int cell = row[e].index();
            // A remote (off-rank) neighbor is never a local box cell.
            if (cell != kRemoteCell && boxOfCell[cell] >= 0) {
                faceDropped[face] = true;
                break;
            }
        }
    }

    // Conformity bookkeeping: (parent, axis, side) -> verified once.
    // The parent's boundary in that direction must consist of exactly one
    // level-zero face that is the parent's *full* logical face (its four
    // lattice corners); otherwise the block boundary is faulted/degenerate
    // and the refined mosaic cannot replace it (unsupported, throw).
    Dune::cpgrid::EntityVariableBase<enum face_tag>& tags0 = faceTag0;
    std::map<std::tuple<int,int,int>, int> verifiedNeighbor;
    const auto outsideNeighborOf = [&](int b, int parent, int axis, int side) {
        const auto key = std::make_tuple(parent, axis, side);
        auto it = verifiedNeighbor.find(key);
        if (it != verifiedNeighbor.end()) {
            return it->second;
        }
        const auto fail = [&](const std::string& what) {
            throw std::logic_error("Refinement box '" + requests[b].name
                                   + "': " + what + " at the block boundary (parent cell "
                                   + std::to_string(parent) + "). Not supported yet.");
        };

        // All faces of the parent in the (axis, side) direction.
        std::vector<int> directionFaces;
        const auto row = cellToFace0[EntityRep<0>(parent, true)];
        for (int e = 0; e < row.size(); ++e) {
            const int face = row[e].index();
            const enum face_tag tag = tags0.get(face);
            if (tag != NNC_FACE && axisOf(tag) == axis
                && row[e].orientation() == (side > 0)) {
                directionFaces.push_back(face);
            }
        }

        int neighbor = -1;
        if (directionFaces.size() > 1) {
            fail("a fault-split face");
        }
        else if (directionFaces.size() == 1) {
            const int face = directionFaces[0];
            // The face must span the parent's full logical face: its corner
            // set must equal the parent's four lattice corners on that side.
            const int fixedBit = (side > 0) ? 1 : 0;
            std::set<int> expected;
            for (int corner = 0; corner < 8; ++corner) {
                const int bit = (axis == 0) ? (corner & 1)
                              : (axis == 1) ? ((corner >> 1) & 1)
                                            : ((corner >> 2) & 1);
                if (bit == fixedBit) {
                    expected.insert(cellToPoint0[parent][corner]);
                }
            }
            std::set<int> actual;
            auto facePoints = faceToPoint0[face];
            for (auto pt = facePoints.begin(); pt != facePoints.end(); ++pt) {
                actual.insert(*pt);
            }
            if (actual != expected) {
                fail("a partial (faulted or degenerate) face");
            }
            const auto cells = faceToCell0[EntityRep<1>(face, true)];
            for (int q = 0; q < cells.size(); ++q) {
                const int cell = cells[q].index();
                if (cell != parent && cell != kRemoteCell) {
                    neighbor = cell;
                }
            }
        }
        // No direction faces: domain boundary, inactive, or off-rank neighbor.

        verifiedNeighbor[key] = neighbor;
        return neighbor;
    };

    // ----- Leaf faces ------------------------------------------------------
    std::vector<SourceRef> leafFaces;
    // outside coarse leaf cell for mosaic faces, aligned with leafFaces (-1 otherwise)
    std::vector<int> mosaicOutside;
    std::vector<int> leafIdxOfFace0(numFaces0, -1);
    for (int face = 0; face < numFaces0; ++face) {
        if (!faceDropped[face]) {
            leafIdxOfFace0[face] = static_cast<int>(leafFaces.size());
            leafFaces.push_back({0, face});
            mosaicOutside.push_back(-1);
        }
    }
    // Face-sharing LGRs: a boundary face of box A toward another refined
    // box B pairs with B's matching boundary face. Both have the same
    // merged corner set (corners already merged), so we key by it. The box
    // seen first owns the leaf face; the partner maps onto it and patches
    // in its cell as the second side (the same mechanism as a coarse
    // mosaic neighbor, only the outside cell is refined).
    std::map<std::vector<int>, int> sharedBoundaryFaceLeaf; // corner set -> owner leaf face
    for (int b = 0; b < numBoxes; ++b) {
        BoxData& box = boxes[b];
        auto& faceToCellL = GridStateWriter::faceToCell(*box.level);
        auto& faceTagL = GridStateWriter::faceTag(*box.level);
        auto& faceToPointL = GridStateWriter::faceToPoint(*box.level);
        Dune::cpgrid::EntityVariableBase<enum face_tag>& tagsL = faceTagL;
        const auto& childToParent = GridStateWriter::childToParent(*box.level);
        box.faceToLeaf.assign(faceToCellL.size(), -1);

        const auto mergedCornerSet = [&](int face) {
            std::vector<int> set;
            auto pts = faceToPointL[face];
            for (auto p = pts.begin(); p != pts.end(); ++p) {
                set.push_back(box.cornerToLeaf[*p]);
            }
            std::sort(set.begin(), set.end());
            return set;
        };

        for (int face = 0; face < faceToCellL.size(); ++face) {
            const auto row = faceToCellL[EntityRep<1>(face, true)];
            int outside = -1;
            if (row.size() == 1) {
                // Boundary face of the refined block: internal hole, domain
                // boundary, block boundary toward a coarse neighbor, or
                // toward another refined box (face-sharing).
                const int cell = row[0].index();
                const bool normalOut = row[0].orientation();
                const int axis = axisOf(tagsL.get(face));
                const int side = normalOut ? 1 : -1;
                const int cart = box.level->globalCell()[cell];
                std::array<int,3> lattice = { cart % box.refinedDims[0],
                                              (cart / box.refinedDims[0]) % box.refinedDims[1],
                                              cart / (box.refinedDims[0]*box.refinedDims[1]) };
                lattice[axis] += side;
                if (lattice[axis] < 0 || lattice[axis] >= box.refinedDims[axis]) {
                    // Block boundary: find the neighbor parent (or none).
                    const int parent = childToParent[cell][1];
                    const int neighbor = outsideNeighborOf(b, parent, axis, side);
                    if (neighbor >= 0 && boxOfCell[neighbor] < 0) {
                        // Toward a coarse (unrefined) neighbor: mosaic.
                        outside = leafIdxOfCell0[neighbor];
                    }
                    else if (neighbor >= 0) {
                        // Toward another refined box: pair the two faces.
                        const auto key = mergedCornerSet(face);
                        const auto it = sharedBoundaryFaceLeaf.find(key);
                        const int thisCellLeaf = leafIdxOfLevelCell[b][cell];
                        if (it != sharedBoundaryFaceLeaf.end()) {
                            // Partner: map onto the owner's face, supply our
                            // cell as its second side, do not emit a face.
                            box.faceToLeaf[face] = it->second;
                            mosaicOutside[it->second] = thisCellLeaf;
                            continue;
                        }
                        // Owner: emit now, partner patches its cell later.
                        sharedBoundaryFaceLeaf.emplace(key, static_cast<int>(leafFaces.size()));
                    }
                }
                // else: internal hole boundary -> stays a boundary face.
            }
            box.faceToLeaf[face] = static_cast<int>(leafFaces.size());
            leafFaces.push_back({b + 1, face});
            mosaicOutside.push_back(outside);
        }
    }
    const int numLeafFaces = static_cast<int>(leafFaces.size());

    // ----- Leaf topology ---------------------------------------------------
    auto& leafFaceToCell = GridStateWriter::faceToCell(*leaf);
    auto& leafFaceToPoint = GridStateWriter::faceToPoint(*leaf);
    auto& leafCellToPoint = GridStateWriter::cellToPoint(*leaf);
    auto& leafCellToFace = GridStateWriter::cellToFace(*leaf);
    auto& leafFaceTag = GridStateWriter::faceTag(*leaf);
    auto& leafFaceNormals = GridStateWriter::faceNormals(*leaf);

    Dune::cpgrid::EntityVariableBase<enum face_tag>& leafTags = leafFaceTag;
    Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& leafNormals = leafFaceNormals;
    leafTags.resize(numLeafFaces);
    leafNormals.resize(numLeafFaces);

    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& leafFaceGeom =
        *(leafGeometry.geomVector(std::integral_constant<int,1>()));
    leafFaceGeom.resize(numLeafFaces);

    std::vector<EntityRep<0>> rowBuffer;
    std::vector<int> pointBuffer;
    for (int face = 0; face < numLeafFaces; ++face) {
        const SourceRef src = leafFaces[face];
        CpGridData& srcGrid = (src.grid == 0) ? level0 : *boxes[src.grid - 1].level;
        auto& srcFaceToCell = (src.grid == 0) ? faceToCell0 : GridStateWriter::faceToCell(srcGrid);
        auto& srcFaceToPoint = (src.grid == 0) ? faceToPoint0 : GridStateWriter::faceToPoint(srcGrid);
        const Dune::cpgrid::EntityVariableBase<enum face_tag>& srcTags =
            (src.grid == 0) ? static_cast<const Dune::cpgrid::EntityVariableBase<enum face_tag>&>(faceTag0)
                            : static_cast<const Dune::cpgrid::EntityVariableBase<enum face_tag>&>(GridStateWriter::faceTag(srcGrid));
        const Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& srcNormals =
            (src.grid == 0) ? static_cast<const Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>&>(faceNormals0)
                            : static_cast<const Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>&>(GridStateWriter::faceNormals(srcGrid));
        const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& srcFaceGeom =
            *(GridStateWriter::geometry(srcGrid).geomVector(std::integral_constant<int,1>()));

        // face_to_cell row (cells mapped to leaf indices; mosaic outside added)
        rowBuffer.clear();
        const auto srcRow = srcFaceToCell[EntityRep<1>(src.index, true)];
        for (int e = 0; e < srcRow.size(); ++e) {
            const int srcCell = srcRow[e].index();
            // Preserve the off-rank sentinel: a face to a remote cell stays
            // a remote-neighbor face on the distributed leaf.
            const int leafCell = (srcCell == kRemoteCell)
                ? kRemoteCell
                : (src.grid == 0)
                    ? leafIdxOfCell0[srcCell]
                    : leafIdxOfLevelCell[src.grid - 1][srcCell];
            rowBuffer.push_back(EntityRep<0>(leafCell, srcRow[e].orientation()));
        }
        if (mosaicOutside[face] >= 0) {
            // The refined cell entry already exists; add the coarse neighbor
            // with opposite orientation, keeping normal-side first.
            const bool refinedNormalOut = srcRow[0].orientation();
            const EntityRep<0> outsideEntry(mosaicOutside[face], !refinedNormalOut);
            if (refinedNormalOut) {
                rowBuffer.push_back(outsideEntry);
            }
            else {
                rowBuffer.insert(rowBuffer.begin(), outsideEntry);
            }
        }
        leafFaceToCell.appendRow(rowBuffer.begin(), rowBuffer.end());

        // face_to_point row (pool indices)
        pointBuffer.clear();
        auto srcPoints = srcFaceToPoint[src.index];
        for (int e = 0; e < srcPoints.size(); ++e) {
            pointBuffer.push_back((src.grid == 0) ? srcPoints[e]
                                                  : boxes[src.grid - 1].cornerToLeaf[srcPoints[e]]);
        }
        leafFaceToPoint.appendRow(pointBuffer.begin(), pointBuffer.end());

        leafTags[face] = srcTags[src.index];
        leafNormals[face] = srcNormals[src.index];
        leafFaceGeom[face] = srcFaceGeom[src.index];
    }

    // cell_to_face_ as the inverse of face_to_cell. Built manually rather
    // than via makeInverseRelation because face_to_cell may carry the
    // off-rank sentinel (distributed grids), which the generic inverse — it
    // sizes the table from the maximum cell index — cannot handle. The
    // orientation convention matches makeInverseRelation: a cell keeps the
    // face's orientation when it is on the face's normal side, else flips it.
    {
        std::vector<std::vector<EntityRep<1>>> cellFaces(numLeafCells);
        for (int face = 0; face < numLeafFaces; ++face) {
            const auto row = leafFaceToCell[EntityRep<1>(face, true)];
            for (int e = 0; e < row.size(); ++e) {
                const int cell = row[e].index();
                if (cell == kRemoteCell) {
                    continue; // remote neighbor: no local cell owns this side
                }
                const EntityRep<1> faceEnt(face, true);
                cellFaces[cell].push_back(row[e].orientation() ? faceEnt : faceEnt.opposite());
            }
        }
        for (int cell = 0; cell < numLeafCells; ++cell) {
            leafCellToFace.appendRow(cellFaces[cell].begin(), cellFaces[cell].end());
        }
    }

    // ----- Leaf cell data ----------------------------------------------------
    leafCellToPoint.resize(numLeafCells);
    std::vector<int> leafGlobalCell(numLeafCells);
    std::vector<std::array<int,2>> leafChildToParent(numLeafCells, std::array<int,2>{-1, -1});
    std::vector<int> leafIdxInParent(numLeafCells, -1);
    std::vector<std::array<int,2>> leafToLevel(numLeafCells);

    for (int cell = 0; cell < numLeafCells; ++cell) {
        const SourceRef src = leafCells[cell];
        if (src.grid == 0) {
            leafCellToPoint[cell] = cellToPoint0[src.index];
            leafGlobalCell[cell] = level0.globalCell()[src.index];
            leafToLevel[cell] = {0, src.index};
        }
        else {
            const BoxData& box = boxes[src.grid - 1];
            const auto& srcPoints = GridStateWriter::cellToPoint(*box.level)[src.index];
            for (int corner = 0; corner < 8; ++corner) {
                leafCellToPoint[cell][corner] = box.cornerToLeaf[srcPoints[corner]];
            }
            const int parent = GridStateWriter::childToParent(*box.level)[src.index][1];
            leafGlobalCell[cell] = level0.globalCell()[parent];
            leafChildToParent[cell] = {0, parent};
            leafIdxInParent[cell] = GridStateWriter::idxInParent(*box.level)[src.index];
            leafToLevel[cell] = {src.grid, src.index};
        }
    }

    // Cell geometries: centers/volumes from the source grids, corner list
    // pointing into the leaf pool via the leaf's cell_to_point_.
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& leafCellGeom =
        *(leafGeometry.geomVector(std::integral_constant<int,0>()));
    leafCellGeom.resize(numLeafCells);
    auto leafCornersPtr = leafGeometry.geomVector(std::integral_constant<int,3>());
    for (int cell = 0; cell < numLeafCells; ++cell) {
        const SourceRef src = leafCells[cell];
        CpGridData& srcGrid = (src.grid == 0) ? level0 : *boxes[src.grid - 1].level;
        const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& srcCellGeom =
            *(GridStateWriter::geometry(srcGrid).geomVector(std::integral_constant<int,0>()));
        const auto& geom = srcCellGeom[src.index];
        leafCellGeom[cell] = Dune::cpgrid::Geometry<3,3>(geom.center(), geom.volume(),
                                                         leafCornersPtr,
                                                         leafCellToPoint[cell].data());
    }

    // ----- Parallel leaf structures (distributed grids only) -------------------
    // Build the leaf's parallel cell index set, remote indices, partition types,
    // and communication interfaces so that overlap leaf cells take part in
    // collective communication (e.g. opm-models' dofTotalVolume border sync,
    // which sums interior volumes and then communicates them to the overlap).
    // Coarse leaf cells inherit their level-zero global index and attribute;
    // refined cells are rank-interior (owner only) and get fresh, globally
    // unique indices in a per-rank range above every coarse global id, so the
    // remote-index matching never pairs them with a cell on another rank.
#if HAVE_MPI
    Dune::Communication<Dune::MPIHelper::MPICommunicator> cc(comm);
    if (cc.size() > 1) {
        using ParallelIndexSet = CpGridData::ParallelIndexSet;
        using AttributeSet = CpGridData::AttributeSet;

        std::vector<int> leafCellGlobalId(numLeafCells, -1);
        std::vector<AttributeSet> leafCellAttr(numLeafCells, AttributeSet::owner);

        int localMaxCoarseGlobal = 0;
        for (const auto& entry : level0.cellIndexSet()) {
            const int l0local = entry.local().local();
            localMaxCoarseGlobal = std::max(localMaxCoarseGlobal, static_cast<int>(entry.global()));
            const int leafC = leafIdxOfCell0[l0local];   // -1 if refined away
            if (leafC >= 0) {
                leafCellGlobalId[leafC] = entry.global();
                leafCellAttr[leafC] = entry.local().attribute();
            }
        }

        // Fresh global ids for refined (interior-only) cells in a per-rank
        // disjoint range above every coarse global id on any rank.
        int numRefinedLocal = 0;
        for (int cell = 0; cell < numLeafCells; ++cell) {
            if (leafCells[cell].grid != 0) ++numRefinedLocal;
        }
        const int globalMaxCoarse = cc.max(localMaxCoarseGlobal);
        std::vector<int> refinedCounts(cc.size(), 0);
        cc.allgather(&numRefinedLocal, 1, refinedCounts.data());
        int refinedNext = globalMaxCoarse + 1;
        for (int r = 0; r < cc.rank(); ++r) {
            refinedNext += refinedCounts[r];
        }
        for (int cell = 0; cell < numLeafCells; ++cell) {
            if (leafCells[cell].grid != 0) {
                leafCellGlobalId[cell] = refinedNext++;
                leafCellAttr[cell] = AttributeSet::owner;
            }
        }

        auto& cis = leaf->cellIndexSet();
        cis.beginResize();
        for (int cell = 0; cell < numLeafCells; ++cell) {
            cis.add(leafCellGlobalId[cell],
                    ParallelIndexSet::LocalIndex(cell, leafCellAttr[cell], true));
        }
        cis.endResize();

        leaf->cellRemoteIndices().template rebuild<false>();
        leaf->computeCellPartitionType();
        leaf->computePointPartitionType();
        leaf->computeCommunicationInterfaces(numLeafCorners);
    }
#endif

    // ----- Leaf metadata -------------------------------------------------------
    GridStateWriter::setLogicalCartesianSize(*leaf, dims0);
    GridStateWriter::setGlobalCell(*leaf, std::move(leafGlobalCell));
    GridStateWriter::setIndexSet(*leaf, numLeafCells, numLeafCorners);
    GridStateWriter::setParentRelations(*leaf, std::move(leafChildToParent), std::move(leafIdxInParent));
    GridStateWriter::setLeafToLevel(*leaf, std::move(leafToLevel));
    GridStateWriter::setCornerHistory(*leaf, std::move(leafCornerHistory));
    GridStateWriter::setRefinementMaxLevel(*leaf, numBoxes);
    for (int b = 0; b < numBoxes; ++b) {
        GridStateWriter::setRefinementMaxLevel(*boxes[b].level, numBoxes);
    }

    // Level zero learns about its children.
    GridStateWriter::setParentToChildren(level0, std::move(parentToChildren));

    return leaf;
}

} // namespace Refinement
} // namespace Opm
