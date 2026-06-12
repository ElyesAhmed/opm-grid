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
    // Boxes must be pairwise separated by at least one cell: touching boxes
    // would require identifying refined corners between two refined level
    // grids, which is not implemented yet.
    for (std::size_t i = 0; i < requests.size(); ++i) {
        for (std::size_t j = i + 1; j < requests.size(); ++j) {
            bool separated = false;
            for (int c = 0; c < 3; ++c) {
                separated = separated
                    || (requests[i].endIJK[c] + 1 <= requests[j].startIJK[c])
                    || (requests[j].endIJK[c] + 1 <= requests[i].startIJK[c]);
            }
            if (!separated) {
                throw std::logic_error("Refinement boxes '" + requests[i].name + "' and '"
                                       + requests[j].name
                                       + "' touch each other. Touching boxes are not supported yet.");
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
