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
#ifndef OPM_GRID_REFINEMENT_GRIDSTATEWRITER_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_GRIDSTATEWRITER_HEADER_INCLUDED

#include <array>
#include <vector>

namespace Dune
{
namespace cpgrid
{
class CpGridData;
}
}

namespace Opm
{
namespace Refinement
{

/// The single write-access point for multilevel/refinement state inside
/// CpGridData. Befriended by CpGridData so the refinement builder can stay
/// outside the grid classes without widening the friend surface: the
/// builder is the only writer of this state (docs/DESIGN-builder.md §5),
/// and every mutation goes through these named operations.
struct GridStateWriter
{
    /// Set the level index a grid occupies in the hierarchy.
    static void setLevel(Dune::cpgrid::CpGridData& grid, int level);

    /// Set the per-parent subdivision factors of a refined level grid.
    static void setCellsPerDim(Dune::cpgrid::CpGridData& grid,
                               const std::array<int,3>& cellsPerDim);

    /// Set the parent relations of a refined level grid:
    /// childToParent[c] = {parent level, parent cell index}, and
    /// idxInParent[c] = lattice index i + j*rx + k*rx*ry within the parent
    /// (the convention geometryInFather() expects).
    static void setParentRelations(Dune::cpgrid::CpGridData& grid,
                                   std::vector<std::array<int,2>> childToParent,
                                   std::vector<int> idxInParent);
};

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_REFINEMENT_GRIDSTATEWRITER_HEADER_INCLUDED
