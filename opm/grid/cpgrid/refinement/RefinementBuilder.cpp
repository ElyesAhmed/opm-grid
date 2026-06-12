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

#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{

std::unique_ptr<Opm::Refinement::Builder>& theBuilder()
{
    static std::unique_ptr<Opm::Refinement::Builder> instance{};
    return instance;
}

bool boxesOverlap(const Opm::Refinement::BlockRefinement& a,
                  const Opm::Refinement::BlockRefinement& b)
{
    for (int c = 0; c < 3; ++c) {
        if (a.endIJK[c] <= b.startIJK[c] || b.endIJK[c] <= a.startIJK[c]) {
            return false;
        }
    }
    return true;
}

} // anonymous namespace

namespace Opm
{
namespace Refinement
{

Builder* builder()
{
    return theBuilder().get();
}

std::unique_ptr<Builder> setBuilder(std::unique_ptr<Builder> newBuilder)
{
    std::unique_ptr<Builder> previous = std::move(theBuilder());
    theBuilder() = std::move(newBuilder);
    return previous;
}

void validateBlockRefinements(const std::vector<BlockRefinement>& requests)
{
    for (const auto& req : requests) {
        for (int c = 0; c < 3; ++c) {
            if (req.startIJK[c] < 0 || req.startIJK[c] >= req.endIJK[c]) {
                throw std::invalid_argument("Invalid IJK box in refinement '" + req.name
                                            + "': end I/J/K must be larger than start I/J/K (and start non-negative).");
            }
            if (req.cellsPerDim[c] < 1) {
                throw std::invalid_argument("Invalid subdivisions in refinement '" + req.name
                                            + "': NX/NY/NZ must be positive.");
            }
        }
    }
    // Pairwise disjointness of boxes refining the same parent grid. Note:
    // upstream never validated this and overlapping CARFIN boxes corrupt
    // the result silently (review Part III.3) - here it is a hard error.
    for (std::size_t i = 0; i < requests.size(); ++i) {
        for (std::size_t j = i + 1; j < requests.size(); ++j) {
            if (requests[i].parentGridName != requests[j].parentGridName) {
                continue;
            }
            if (boxesOverlap(requests[i], requests[j])) {
                throw std::invalid_argument("Refinement boxes '" + requests[i].name + "' and '"
                                            + requests[j].name + "' overlap. CARFIN boxes must be disjoint.");
            }
        }
    }
}

} // namespace Refinement
} // namespace Opm
