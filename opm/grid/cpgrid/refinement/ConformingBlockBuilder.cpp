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

#include <opm/grid/cpgrid/refinement/ConformingBlockBuilder.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/refinement/EdgeConformal.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>
#include <opm/grid/cpgrid/refinement/LeafGridAssembler.hpp>
#include <opm/grid/cpgrid/refinement/LevelGridAssembler.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <stdexcept>
#include <string>
#include <utility>

namespace Opm
{
namespace Refinement
{

ConformingBlockBuilder::ConformingBlockBuilder(const std::array<int,3>& parentDims,
                                               std::vector<double> coord,
                                               std::vector<double> zcorn,
                                               std::vector<int> actnum,
                                               bool edgeConformal)
    : dims_(parentDims)
    , coord_(std::move(coord))
    , zcorn_(std::move(zcorn))
    , actnum_(std::move(actnum))
    , edgeConformal_(edgeConformal)
{
}

namespace
{

// Per-box presence on this rank (distributed runs). The rank-interior model
// (PLAN Track 1 step 6): a box is either fully owned here (all its cells
// local and interior) or absent (no cell local, not even as overlap).
// Anything in between means the box crosses a rank boundary / touches the
// overlap, which is not supported.
enum class BoxPresence { Owned, Absent };

BoxPresence classifyBox(const Dune::cpgrid::CpGridData& level0,
                        const std::array<int,3>& dims,
                        const Opm::Refinement::BlockRefinement& req)
{
    const int boxCells = (req.endIJK[0] - req.startIJK[0])
                       * (req.endIJK[1] - req.startIJK[1])
                       * (req.endIJK[2] - req.startIJK[2]);
    int interior = 0;
    int nonInterior = 0;
    const auto& globalCell = level0.globalCell();
    for (int c = 0; c < level0.size(0); ++c) {
        const int cart = globalCell[c];
        const int i = cart % dims[0];
        const int j = (cart / dims[0]) % dims[1];
        const int k = cart / (dims[0]*dims[1]);
        if (i >= req.startIJK[0] && i < req.endIJK[0]
            && j >= req.startIJK[1] && j < req.endIJK[1]
            && k >= req.startIJK[2] && k < req.endIJK[2]) {
            const Dune::cpgrid::Entity<0> e(level0, c, true);
            if (e.partitionType() == Dune::InteriorEntity) {
                ++interior;
            } else {
                ++nonInterior;
            }
        }
    }
    if (interior == 0 && nonInterior == 0) {
        return BoxPresence::Absent;
    }
    if (interior == boxCells && nonInterior == 0) {
        return BoxPresence::Owned;
    }
    throw std::logic_error("Refinement box '" + req.name + "' is split across MPI ranks or "
                           "touches the overlap region. Only rank-interior LGRs are supported; "
                           "use CpGrid::setPartitionCellGroups to keep the box on one rank.");
}

// refine-before-redistribute classification: before load balancing the global
// grid lives on rank 0 (all cells) and is empty on the other ranks. Classify a
// box purely by cell *presence* (ignoring partition type, which is not yet set):
// the rank that holds all of the box's cells refines it; ranks with none carry
// an empty placeholder level grid. A partial presence (some but not all cells)
// is not expected here and is rejected.
BoxPresence boxCellPresence(const Dune::cpgrid::CpGridData& level0,
                            const std::array<int,3>& dims,
                            const Opm::Refinement::BlockRefinement& req)
{
    const int boxCells = (req.endIJK[0] - req.startIJK[0])
                       * (req.endIJK[1] - req.startIJK[1])
                       * (req.endIJK[2] - req.startIJK[2]);
    int present = 0;
    const auto& globalCell = level0.globalCell();
    for (int c = 0; c < level0.size(0); ++c) {
        const int cart = globalCell[c];
        const int i = cart % dims[0];
        const int j = (cart / dims[0]) % dims[1];
        const int k = cart / (dims[0]*dims[1]);
        if (i >= req.startIJK[0] && i < req.endIJK[0]
            && j >= req.startIJK[1] && j < req.endIJK[1]
            && k >= req.startIJK[2] && k < req.endIJK[2]) {
            ++present;
        }
    }
    if (present == 0) {
        return BoxPresence::Absent;
    }
    if (present == boxCells) {
        return BoxPresence::Owned;
    }
    throw std::logic_error("Refinement box '" + req.name + "' is only partially present on "
                           "this rank before load balancing; refine-before-redistribute "
                           "expects the full global grid on a single rank.");
}

} // anonymous namespace

void ConformingBlockBuilder::build(Dune::CpGrid& grid,
                                   const std::vector<BlockRefinement>& requests)
{
    if (grid.maxLevel() != 0) {
        throw std::logic_error("ConformingBlockBuilder requires an unrefined starting grid.");
    }
    for (const auto& req : requests) {
        if (req.parentGridName != "GLOBAL") {
            throw std::logic_error("Nested refinement ('" + req.name + "' with parent '"
                                   + req.parentGridName + "') is not supported yet.");
        }
    }
    // Conformity check for touching boxes. With disjoint boxes, a pair
    // either is fully separated, shares an edge/corner, or shares a 2D
    // face. The shared refined entities (corners along a shared edge,
    // corners+faces on a shared face) coincide bitwise only if the
    // subdivision factors match in the *shared* directions: the overlap
    // dimension(s). The touch dimension(s) are parent-cell boundaries where
    // corners are identified with level zero regardless. Subdivisions in a
    // dimension where the boxes neither touch nor overlap are unconstrained
    // (e.g. the out-of-face direction of two face-sharing boxes).
    for (std::size_t i = 0; i < requests.size(); ++i) {
        for (std::size_t j = i + 1; j < requests.size(); ++j) {
            const auto& a = requests[i];
            const auto& b = requests[j];
            // The boxes share an entity (corner/edge/face) only if they
            // touch-or-overlap in every dimension; a gap in any dimension
            // means they are fully separated and impose no constraint.
            bool interacts = true;
            for (int c = 0; c < 3; ++c) {
                const bool gap = (a.endIJK[c] < b.startIJK[c]) || (b.endIJK[c] < a.startIJK[c]);
                interacts = interacts && !gap;
            }
            if (!interacts) {
                continue;
            }
            for (int c = 0; c < 3; ++c) {
                const bool overlap = (a.startIJK[c] < b.endIJK[c]) && (b.startIJK[c] < a.endIJK[c]);
                if (overlap && a.cellsPerDim[c] != b.cellsPerDim[c]) {
                    throw std::logic_error("Refinement boxes '" + a.name + "' and '" + b.name
                                           + "' meet with non-matching subdivisions in direction "
                                           + std::to_string(c) + ". Not supported.");
                }
            }
        }
    }

    auto& storage = grid.currentData();
    // The refinement is "distributed" (rank-interior: only the owning rank
    // refines each box) only once the grid has actually been scattered. Before
    // load balancing - the refine-before-redistribute path - the grid is still
    // the replicated global grid; refine it serially on every rank instead, so
    // the builder refines every cell (no rank-interior classifyBox) and the
    // leaf assembly does no collectives, producing the full refined grid a
    // serial run yields, ready to be distributed afterwards.
    const bool distributed = grid.comm().size() > 1 && grid.isDistributed();
    // Capture the grid's actual (world) communicator before any new level or
    // leaf grid is pushed - CpGrid::comm() reads the back() entry of the
    // storage, which is the coarse level-zero grid at this point and therefore
    // still carries the world communicator on every rank.
    const auto worldComm = grid.comm();
    // Use a self-communicator for serial refinement (refine-before-redistribute
    // and the self-comm output/reference grid) so no operation is collective
    // over the world communicator; the distributed simulation grid keeps the
    // world communicator.
    Dune::MPIHelper::MPICommunicator comm = Dune::MPIHelper::getLocalCommunicator();
    if (distributed) {
        comm = grid.comm();
    }
    // refine-before-redistribute: comm().size() > 1 but the grid is not yet
    // scattered. The global grid is on rank 0 and empty on the others, so we
    // classify boxes by cell presence (rank 0 refines, the rest get empty
    // placeholders) and refine serially with the self-communicator above.
    const bool refineBefore = grid.comm().size() > 1 && !grid.isDistributed();

    for (std::size_t b = 0; b < requests.size(); ++b) {
        // In a distributed run, only the rank that owns the (rank-interior)
        // box refines it; other ranks carry an empty placeholder level grid.
        BoxPresence presence = BoxPresence::Owned;
        if (distributed) {
            presence = classifyBox(*storage[0], dims_, requests[b]);
        } else if (refineBefore) {
            presence = boxCellPresence(*storage[0], dims_, requests[b]);
        }
        auto level = (presence == BoxPresence::Owned)
            ? assembleBlockLevelGrid(*storage[0], dims_,
                                     coord_.data(), zcorn_.data(),
                                     actnum_.empty() ? nullptr : actnum_.data(),
                                     requests[b], static_cast<int>(b) + 1,
                                     storage, comm)
            : assembleEmptyLevelGrid(requests[b], static_cast<int>(b) + 1, storage, comm);
        storage.push_back(std::move(level));
    }

    auto leaf = assembleLeafGrid(storage, requests, dims_,
                                 coord_.data(), zcorn_.data(),
                                 actnum_.empty() ? nullptr : actnum_.data(),
                                 comm);
    if (edgeConformal_) {
        edgeConformalizeLeaf(*leaf);
    }
    storage.push_back(std::move(leaf));

    if (refineBefore) {
        // The new level grids (indices 1..N) and the leaf (back) were assembled
        // with a self-communicator so the serial refinement does no collectives.
        // Restore the world communicator on them: CpGrid::comm() reads the leaf,
        // so without this the refined grid would report a size-1 communicator
        // and the subsequent load balancing would treat it as a serial grid and
        // skip distributing the refined leaf. Level zero (index 0) is the
        // pre-existing coarse grid and already carries the world communicator.
        for (std::size_t i = 1; i < storage.size(); ++i) {
            GridStateWriter::setCommunicator(*storage[i], worldComm);
        }
    }
}

} // namespace Refinement
} // namespace Opm
