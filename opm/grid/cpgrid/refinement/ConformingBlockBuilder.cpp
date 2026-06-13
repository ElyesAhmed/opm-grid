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
                                               std::vector<int> actnum)
    : dims_(parentDims)
    , coord_(std::move(coord))
    , zcorn_(std::move(zcorn))
    , actnum_(std::move(actnum))
{
}

void ConformingBlockBuilder::build(Dune::CpGrid& grid,
                                   const std::vector<BlockRefinement>& requests)
{
    if (grid.comm().size() > 1) {
        throw std::logic_error("ConformingBlockBuilder supports serial runs only, yet.");
    }
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
    const auto comm = Dune::MPIHelper::getCommunicator();

    for (std::size_t b = 0; b < requests.size(); ++b) {
        auto level = assembleBlockLevelGrid(*storage[0], dims_,
                                            coord_.data(), zcorn_.data(),
                                            actnum_.empty() ? nullptr : actnum_.data(),
                                            requests[b], static_cast<int>(b) + 1,
                                            storage, comm);
        storage.push_back(std::move(level));
    }

    auto leaf = assembleLeafGrid(storage, requests, comm);
    storage.push_back(std::move(leaf));
}

} // namespace Refinement
} // namespace Opm
