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

#include <opm/grid/cpgrid/refinement/LevelGridAssembler.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/refinement/GrdeclRefinement.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace Opm
{
namespace Refinement
{

std::shared_ptr<Dune::cpgrid::CpGridData>
assembleBlockLevelGrid(const Dune::cpgrid::CpGridData& parentGrid,
                       const std::array<int,3>& parentDims,
                       const double* coord,
                       const double* zcorn,
                       const int* actnum,
                       const BlockRefinement& request,
                       int levelIndex,
                       int parentLevel,
                       std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& levelStorage,
                       Dune::MPIHelper::MPICommunicator comm,
                       RefinedBlockGrdecl* outRefined)
{
    // Stages 2+3: resample the block geometry onto sub-pillars.
    const RefinedBlockGrdecl refined = refineBlock(parentDims, coord, zcorn, actnum, request);

    // Stage 4: process the refined description with the corner-point
    // preprocessor. Faults/pinch-outs inside the block are matched here.
    // Refined level grids are rank-local (rank-interior LGR model): use a
    // self-communicator so none of their operations are collective over the
    // distributed communicator (avoids deadlock when only the owning rank
    // refines). The distributed leaf keeps the real communicator.
    (void)comm;
    auto level = std::make_shared<Dune::cpgrid::CpGridData>(
        Dune::MPIHelper::getLocalCommunicator(), levelStorage);
    grdecl raw;
    raw.dims[0] = refined.dims[0];
    raw.dims[1] = refined.dims[1];
    raw.dims[2] = refined.dims[2];
    raw.coord = refined.coord.data();
    raw.zcorn = refined.zcorn.data();
    raw.actnum = refined.actnum.data();

    std::array<std::set<std::pair<int,int>>, 2> nnc;
    level->processEclipseFormat(raw,
#if HAVE_OPM_COMMON
                                nullptr,
#endif
                                nnc,
                                /* remove_ij_boundary = */ false,
                                /* turn_normals = */ false,
                                /* pinchActive = */ false,
                                /* tolerance_unique_points = */ 0.0,
                                /* edge_conformal = */ false);

    // Stage 6a: parent relations. The processed grid's global_cell_ holds
    // the refined-local Cartesian index (the CARFIN convention for level
    // grids), which encodes the parent and the position inside it.
    const auto& [rx, ry, rz] = request.cellsPerDim;
    const auto& [nx, ny, nz] = parentDims;

    // Invert the parent grid's cartesian -> compressed mapping. For a
    // top-level box the parent is level zero; for a nested box it is the LGR
    // level grid being refined, whose globalCell() holds the refined-local
    // Cartesian index (the same CARFIN convention), so the same inversion
    // works in the parent's local index space.
    std::vector<int> parentCompressed(static_cast<std::size_t>(nx)*ny*nz, -1);
    const auto& parentGlobalCell = parentGrid.globalCell();
    for (std::size_t c = 0; c < parentGlobalCell.size(); ++c) {
        parentCompressed[parentGlobalCell[c]] = static_cast<int>(c);
    }

    const int numRefinedCells = level->size(0);
    std::vector<std::array<int,2>> childToParent(numRefinedCells);
    std::vector<int> idxInParent(numRefinedCells);

    for (int cell = 0; cell < numRefinedCells; ++cell) {
        const int refinedCart = level->globalCell()[cell];
        const int ir = refinedCart % refined.dims[0];
        const int jr = (refinedCart / refined.dims[0]) % refined.dims[1];
        const int kr = refinedCart / (refined.dims[0]*refined.dims[1]);

        const int ci = request.startIJK[0] + ir / rx;
        const int cj = request.startIJK[1] + jr / ry;
        const int ck = request.startIJK[2] + kr / rz;
        const int parentCart = ci + nx*cj + nx*ny*ck;
        const int parentIdx = parentCompressed[parentCart];
        if (parentIdx < 0) {
            // Children inherit ACTNUM, so an active child of an inactive
            // parent indicates an internal inconsistency.
            throw std::logic_error("Refined cell in '" + request.name
                                   + "' has no active parent cell.");
        }
        childToParent[cell] = {parentLevel, parentIdx};
        idxInParent[cell] = (ir % rx) + (jr % ry)*rx + (kr % rz)*rx*ry;
    }

    GridStateWriter::setLevel(*level, levelIndex);
    GridStateWriter::setCellsPerDim(*level, request.cellsPerDim);
    GridStateWriter::setParentRelations(*level, std::move(childToParent), std::move(idxInParent));

    // Expose the resampled description so a nested child can be refined from
    // this level's geometry (its parent) rather than the global grid.
    if (outRefined != nullptr) {
        *outRefined = refined;
    }

    return level;
}

std::shared_ptr<Dune::cpgrid::CpGridData>
assembleEmptyLevelGrid(const BlockRefinement& request,
                       int levelIndex,
                       std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& levelStorage,
                       Dune::MPIHelper::MPICommunicator comm)
{
    (void)comm;
    auto level = std::make_shared<Dune::cpgrid::CpGridData>(
        Dune::MPIHelper::getLocalCommunicator(), levelStorage);
    const std::array<int,3> refinedDims = {
        (request.endIJK[0] - request.startIJK[0]) * request.cellsPerDim[0],
        (request.endIJK[1] - request.startIJK[1]) * request.cellsPerDim[1],
        (request.endIJK[2] - request.startIJK[2]) * request.cellsPerDim[2] };
    GridStateWriter::setLevel(*level, levelIndex);
    GridStateWriter::setCellsPerDim(*level, request.cellsPerDim);
    GridStateWriter::setLogicalCartesianSize(*level, refinedDims);
    GridStateWriter::setGlobalCell(*level, {});
    GridStateWriter::setIndexSet(*level, 0, 0);
    GridStateWriter::setParentRelations(*level, {}, {});
    return level;
}

} // namespace Refinement
} // namespace Opm
