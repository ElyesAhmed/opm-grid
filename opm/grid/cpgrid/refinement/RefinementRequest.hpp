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
#ifndef OPM_GRID_REFINEMENT_REQUEST_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_REQUEST_HEADER_INCLUDED

#include <array>
#include <string>
#include <vector>

namespace Opm
{
namespace Refinement
{

/// One block-shaped static refinement request (one CARFIN box).
///
/// IJK indices refer to the logical Cartesian space of the parent grid
/// ("GLOBAL" for level zero). The box is half-open in cell terms: the last
/// refined cell column is endIJK - 1, matching CpGrid::addLgrsUpdateLeafView.
struct BlockRefinement
{
    std::string name;
    std::string parentGridName{"GLOBAL"};
    std::array<int,3> cellsPerDim{};
    std::array<int,3> startIJK{};
    std::array<int,3> endIJK{};
};

/// Validate a set of block requests: matching sizes are the caller's
/// responsibility; this checks per-box consistency (start < end, positive
/// subdivisions) and pairwise disjointness of boxes sharing a parent grid.
/// Throws std::invalid_argument with a message naming the offending box(es).
void validateBlockRefinements(const std::vector<BlockRefinement>& requests);

/// Result of coalesceBlockRefinements(): the merged request list plus, for
/// every input box, the index of the merged box that now covers it.
struct CoalescedBlockRefinements
{
    std::vector<BlockRefinement> merged;
    std::vector<int> originalToMerged;   // size == input.size()
};

/// Merge axis-abutting refinement boxes that share a parent grid AND have
/// identical subdivision factors into single boxes whose union is a rectangular
/// block. Each merged box becomes ONE LGR level instead of one level per input
/// box -- N mutually adjacent 1-cell CARFIN boxes otherwise stack into N nested
/// levels, whose coarse<->fine leaf transmissibilities are badly inconsistent
/// (review 2026-09-10). The input must already be valid (disjoint). The merged
/// name is the first contributing box's name.
CoalescedBlockRefinements
coalesceBlockRefinements(const std::vector<BlockRefinement>& requests);

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_REFINEMENT_REQUEST_HEADER_INCLUDED
