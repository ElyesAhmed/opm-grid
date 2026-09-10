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

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

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

namespace
{

// Two disjoint boxes with the same parent+factors merge into a rectangular
// block iff they are equal on two axes and abut (or touch) on the third.
bool mergeableAlongAxis(const BlockRefinement& a, const BlockRefinement& b, int axis)
{
    if (a.parentGridName != b.parentGridName || a.cellsPerDim != b.cellsPerDim) {
        return false;
    }
    for (int c = 0; c < 3; ++c) {
        if (c == axis) {
            continue;
        }
        if (a.startIJK[c] != b.startIJK[c] || a.endIJK[c] != b.endIJK[c]) {
            return false;
        }
    }
    return a.endIJK[axis] >= b.startIJK[axis] && b.endIJK[axis] >= a.startIJK[axis];
}

} // anonymous namespace

CoalescedBlockRefinements
coalesceBlockRefinements(const std::vector<BlockRefinement>& requests)
{
    CoalescedBlockRefinements out;
    out.merged = requests;
    out.originalToMerged.resize(requests.size());
    for (std::size_t i = 0; i < requests.size(); ++i) {
        out.originalToMerged[i] = static_cast<int>(i);
    }

    // Track, for every current merged box, which input boxes it now covers.
    std::vector<std::vector<int>> members(requests.size());
    for (std::size_t i = 0; i < requests.size(); ++i) {
        members[i].push_back(static_cast<int>(i));
    }

    // Axis sweeps to a fixed point. Each sweep sorts so that boxes mergeable
    // along `axis` land next to each other, then a single linear pass absorbs
    // consecutive abutting boxes. O(n log n) per sweep; the list only shrinks.
    bool changed = true;
    while (changed) {
        changed = false;
        for (int axis = 0; axis < 3; ++axis) {
            const int u = (axis + 1) % 3;
            const int v = (axis + 2) % 3;
            std::vector<std::size_t> order(out.merged.size());
            for (std::size_t i = 0; i < order.size(); ++i) {
                order[i] = i;
            }
            std::sort(order.begin(), order.end(),
                [&](std::size_t x, std::size_t y) {
                    const auto& a = out.merged[x];
                    const auto& b = out.merged[y];
                    auto ka = std::make_tuple(a.parentGridName, a.cellsPerDim[0], a.cellsPerDim[1],
                                              a.cellsPerDim[2], a.startIJK[u], a.endIJK[u],
                                              a.startIJK[v], a.endIJK[v], a.startIJK[axis]);
                    auto kb = std::make_tuple(b.parentGridName, b.cellsPerDim[0], b.cellsPerDim[1],
                                              b.cellsPerDim[2], b.startIJK[u], b.endIJK[u],
                                              b.startIJK[v], b.endIJK[v], b.startIJK[axis]);
                    return ka < kb;
                });

            std::vector<BlockRefinement> next;
            std::vector<std::vector<int>> nextMembers;
            next.reserve(out.merged.size());
            for (std::size_t oi = 0; oi < order.size(); ++oi) {
                const std::size_t idx = order[oi];
                if (!next.empty() && mergeableAlongAxis(next.back(), out.merged[idx], axis)) {
                    next.back().startIJK[axis] =
                        std::min(next.back().startIJK[axis], out.merged[idx].startIJK[axis]);
                    next.back().endIJK[axis] =
                        std::max(next.back().endIJK[axis], out.merged[idx].endIJK[axis]);
                    nextMembers.back().insert(nextMembers.back().end(),
                                              members[idx].begin(), members[idx].end());
                    changed = true;
                }
                else {
                    next.push_back(out.merged[idx]);
                    nextMembers.push_back(members[idx]);
                }
            }
            out.merged = std::move(next);
            members = std::move(nextMembers);
        }
    }

    // Restore an input-consistent order: a merged box sits where its
    // lowest-indexed contributing input box was. This keeps any
    // parent-before-child ordering the caller relied on.
    std::vector<std::size_t> perm(members.size());
    for (std::size_t i = 0; i < perm.size(); ++i) {
        perm[i] = i;
    }
    std::sort(perm.begin(), perm.end(), [&](std::size_t x, std::size_t y) {
        return *std::min_element(members[x].begin(), members[x].end())
             < *std::min_element(members[y].begin(), members[y].end());
    });
    std::vector<BlockRefinement> orderedMerged;
    orderedMerged.reserve(perm.size());
    for (std::size_t i = 0; i < perm.size(); ++i) {
        orderedMerged.push_back(out.merged[perm[i]]);
        for (int orig : members[perm[i]]) {
            out.originalToMerged[orig] = static_cast<int>(i);
        }
    }
    out.merged = std::move(orderedMerged);
    return out;
}

} // namespace Refinement
} // namespace Opm
