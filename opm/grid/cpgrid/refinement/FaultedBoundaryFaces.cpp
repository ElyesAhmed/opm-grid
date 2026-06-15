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

#include <opm/grid/cpgrid/refinement/FaultedBoundaryFaces.hpp>

#include <opm/grid/cpgrid/refinement/GrdeclRefinement.hpp>
#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <array>

namespace Opm
{
namespace Refinement
{

std::vector<BoundaryConnection>
faultedBoundaryConnections(const std::array<int,3>& parentDims,
                           const double* coord,
                           const double* zcorn,
                           const int* actnum,
                           const std::array<int,3>& boxStartIJK,
                           const std::array<int,3>& boxEndIJK,
                           const std::array<int,3>& cellsPerDim,
                           int axis,
                           int side,
                           bool edgeConformal)
{
    std::vector<BoundaryConnection> result;

    // Extend the box by one parent cell on the (axis, side) side: the coarse
    // "shell". Refining the box together with this shell and processing the
    // result reproduces the box's boundary faces against the shell, split at
    // any fault by the corner-point processor.
    BlockRefinement req;
    req.name = "FAULTBND";
    req.cellsPerDim = cellsPerDim;
    req.startIJK = boxStartIJK;
    req.endIJK = boxEndIJK;
    if (side < 0) {
        req.startIJK[axis] -= 1;
        if (req.startIJK[axis] < 0) {
            return result; // shell off-grid: domain boundary, no connections
        }
    }
    else {
        req.endIJK[axis] += 1;
        if (req.endIJK[axis] > parentDims[axis]) {
            return result;
        }
    }

    const RefinedBlockGrdecl refined = refineBlock(parentDims, coord, zcorn, actnum, req);
    const std::array<int,3> mdims = refined.dims;

    grdecl raw;
    raw.dims[0] = mdims[0]; raw.dims[1] = mdims[1]; raw.dims[2] = mdims[2];
    raw.coord = refined.coord.data();
    raw.zcorn = refined.zcorn.data();
    raw.actnum = refined.actnum.empty() ? nullptr : refined.actnum.data();

    struct processed_grid out;
    if (!process_grdecl(/*pinchActive=*/1, edgeConformal ? 1 : 0, 1e-6, &raw, nullptr, &out)) {
        return result;
    }

    // The shell occupies the first cellsPerDim[axis] cells along `axis` for
    // side<0, the last cellsPerDim[axis] for side>0; the box is the rest.
    const int f = cellsPerDim[axis];
    const int shellLo = (side < 0) ? 0 : (mdims[axis] - f);
    const int shellHi = shellLo + f;            // [shellLo, shellHi) are shell cells
    const enum face_tag axisFaceTag = (axis == 0) ? I_FACE
                                    : (axis == 1) ? J_FACE
                                                  : K_FACE;

    const auto latticeOf = [&](int activeCell) {
        const int g = out.local_cell_index[activeCell];   // active -> logical (mini) index
        return std::array<int,3>{ g % mdims[0],
                                  (g / mdims[0]) % mdims[1],
                                  g / (mdims[0] * mdims[1]) };
    };

    for (unsigned face = 0; face < out.number_of_faces; ++face) {
        if (out.face_tag[face] != axisFaceTag) {
            continue;
        }
        const int a = out.face_neighbors[2*face];
        const int b = out.face_neighbors[2*face + 1];
        if (a < 0 || b < 0) {
            continue; // domain boundary on this face
        }
        const std::array<int,3> la = latticeOf(a);
        const std::array<int,3> lb = latticeOf(b);
        const bool aShell = (la[axis] >= shellLo && la[axis] < shellHi);
        const bool bShell = (lb[axis] >= shellLo && lb[axis] < shellHi);
        if (aShell == bShell) {
            continue; // both shell or both box: not a box-boundary connection
        }
        const std::array<int,3> shellL = aShell ? la : lb;
        const std::array<int,3> boxL   = aShell ? lb : la;

        // Box cell lattice in the box's own refined frame: for side<0 the box
        // is offset by f along axis; for side>0 it starts at 0.
        std::array<int,3> boxFrame = boxL;
        if (side < 0) {
            boxFrame[axis] -= f;
        }

        // Coarse neighbour parent Cartesian: the shell cell's parent.
        std::array<int,3> parentIJK{};
        for (int d = 0; d < 3; ++d) {
            parentIJK[d] = req.startIJK[d] + shellL[d] / cellsPerDim[d];
        }
        const int coarseCart = parentIJK[0]
                             + parentDims[0] * parentIJK[1]
                             + parentDims[0] * parentDims[1] * parentIJK[2];

        BoundaryConnection bc;
        bc.boxCell = boxFrame;
        bc.coarseNeighborCart = coarseCart;
        for (unsigned n = out.face_node_ptr[face]; n < out.face_node_ptr[face + 1]; ++n) {
            const int node = out.face_nodes[n];
            bc.faceNodes.push_back({ out.node_coordinates[3*node],
                                     out.node_coordinates[3*node + 1],
                                     out.node_coordinates[3*node + 2] });
        }
        result.push_back(std::move(bc));
    }

    free_processed_grid(&out);
    return result;
}

} // namespace Refinement
} // namespace Opm
