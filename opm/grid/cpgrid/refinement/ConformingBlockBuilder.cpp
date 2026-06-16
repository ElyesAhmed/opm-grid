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
    // Use the grid's own communicator (not the global MPIHelper one) so that a
    // grid carrying a self-communicator - e.g. a serial output/reference grid
    // built on the I/O rank of a parallel run - is refined purely locally
    // (no collective calls).  For the distributed simulation grid this is the
    // world communicator, i.e. unchanged behaviour.
    const Dune::MPIHelper::MPICommunicator comm = grid.comm();
    const bool distributed = grid.comm().size() > 1;

    for (std::size_t b = 0; b < requests.size(); ++b) {
        // In a distributed run, only the rank that owns the (rank-interior)
        // box refines it; other ranks carry an empty placeholder level grid.
        const auto presence = distributed
            ? classifyBox(*storage[0], dims_, requests[b])
            : BoxPresence::Owned;
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
}

} // namespace Refinement
} // namespace Opm
