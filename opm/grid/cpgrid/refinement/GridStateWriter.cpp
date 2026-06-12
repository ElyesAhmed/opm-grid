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

#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>

#include <opm/grid/cpgrid/CpGridData.hpp>

#include <utility>

namespace Opm
{
namespace Refinement
{

void GridStateWriter::setLevel(Dune::cpgrid::CpGridData& grid, int level)
{
    grid.level_ = level;
}

void GridStateWriter::setCellsPerDim(Dune::cpgrid::CpGridData& grid,
                                     const std::array<int,3>& cellsPerDim)
{
    grid.cells_per_dim_ = cellsPerDim;
}

void GridStateWriter::setParentRelations(Dune::cpgrid::CpGridData& grid,
                                         std::vector<std::array<int,2>> childToParent,
                                         std::vector<int> idxInParent)
{
    grid.child_to_parent_cells_ = std::move(childToParent);
    grid.cell_to_idxInParentCell_ = std::move(idxInParent);
}

} // namespace Refinement
} // namespace Opm
