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
#ifndef OPM_GRID_REFINEMENT_LEVELGRIDASSEMBLER_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_LEVELGRIDASSEMBLER_HEADER_INCLUDED

#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <memory>
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

/// Stage 4+6a of the build pipeline (docs/DESIGN-builder.md D4): turn one
/// block request into a fully processed refined *level grid*.
///
/// The block is resampled onto sub-pillars (refineBlock), run through the
/// corner-point preprocessor (so faults and pinch-outs inside the block are
/// matched by the same code that handles them on level zero), and equipped
/// with the parent relations the rest of CpGrid reads (father(),
/// geometryInFather(), CARFIN-local Cartesian indices).
///
/// The returned grid has level_ = levelIndex and its level_data_ptr_ bound
/// to levelStorage; the caller is responsible for placing it at position
/// levelIndex in that vector (and level zero at position 0).
///
/// @param level0       The grid the parents live in (currently level zero).
/// @param parentDims   Cartesian dimensions of the parent description.
/// @param coord,zcorn  Parent COORD/ZCORN arrays.
/// @param actnum       Parent ACTNUM or nullptr (all active).
/// @param request      The block to refine.
/// @param levelIndex   Index this level grid will occupy in the hierarchy.
/// @param levelStorage The hierarchy vector the grid will live in.
/// @param comm         Communicator for the new grid object.
std::shared_ptr<Dune::cpgrid::CpGridData>
assembleBlockLevelGrid(const Dune::cpgrid::CpGridData& level0,
                       const std::array<int,3>& parentDims,
                       const double* coord,
                       const double* zcorn,
                       const int* actnum,
                       const BlockRefinement& request,
                       int levelIndex,
                       std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& levelStorage,
                       Dune::MPIHelper::MPICommunicator comm);

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_REFINEMENT_LEVELGRIDASSEMBLER_HEADER_INCLUDED
