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

#include <opm/grid/cpgrid/refinement/GrdeclRefinement.hpp>

#include <stdexcept>
#include <string>

namespace
{

// Lateral position of a refined pillar/corner line within the parent block:
// which parent column it belongs to and the fractional coordinate in it.
// The last refined line is assigned to the last column with fraction 1 so
// that every (column, fraction) pair is well defined.
struct LateralPos
{
    int cell;       // parent cell index in this direction (absolute)
    double frac;    // in [0, 1]
};

LateralPos lateralPos(int refinedIdx, int factor, int start, int numParentCells)
{
    if (refinedIdx == numParentCells * factor) {
        return { start + numParentCells - 1, 1.0 };
    }
    return { start + refinedIdx / factor,
             static_cast<double>(refinedIdx % factor) / factor };
}

} // anonymous namespace

namespace Opm
{
namespace Refinement
{

RefinedBlockGrdecl refineBlock(const std::array<int,3>& parentDims,
                               const double* coord,
                               const double* zcorn,
                               const int* actnum,
                               const BlockRefinement& request)
{
    const auto& [nx, ny, nz] = parentDims;
    for (int c = 0; c < 3; ++c) {
        if (request.startIJK[c] < 0 || request.endIJK[c] > parentDims[c]
            || request.startIJK[c] >= request.endIJK[c]) {
            throw std::invalid_argument("Refinement box '" + request.name
                                        + "' does not fit inside the parent dimensions.");
        }
        if (request.cellsPerDim[c] < 1) {
            throw std::invalid_argument("Refinement '" + request.name
                                        + "' has non-positive subdivisions.");
        }
    }

    const std::array<int,3> boxDims = { request.endIJK[0] - request.startIJK[0],
                                        request.endIJK[1] - request.startIJK[1],
                                        request.endIJK[2] - request.startIJK[2] };
    const std::array<int,3>& factors = request.cellsPerDim;

    RefinedBlockGrdecl out;
    out.dims = { boxDims[0]*factors[0], boxDims[1]*factors[1], boxDims[2]*factors[2] };

    // --- COORD: sub-pillars ---------------------------------------------
    // Endpoint-wise bilinear interpolation of the four parent pillars
    // surrounding the lateral position. On parent pillar positions the
    // weights collapse and the parent pillar is reproduced exactly, so
    // adjacent parent columns see identical sub-pillars.
    const auto parentPillar = [&](int i, int j) {
        return coord + 6*(static_cast<std::size_t>(j)*(nx + 1) + i);
    };

    out.coord.resize(6 * static_cast<std::size_t>(out.dims[0] + 1) * (out.dims[1] + 1));
    for (int jr = 0; jr <= out.dims[1]; ++jr) {
        const auto [cj, b] = lateralPos(jr, factors[1], request.startIJK[1], boxDims[1]);
        for (int ir = 0; ir <= out.dims[0]; ++ir) {
            const auto [ci, a] = lateralPos(ir, factors[0], request.startIJK[0], boxDims[0]);

            const double* p00 = parentPillar(ci,     cj);
            const double* p10 = parentPillar(ci + 1, cj);
            const double* p01 = parentPillar(ci,     cj + 1);
            const double* p11 = parentPillar(ci + 1, cj + 1);

            double* sub = &out.coord[6*(static_cast<std::size_t>(jr)*(out.dims[0] + 1) + ir)];
            for (int comp = 0; comp < 6; ++comp) {
                sub[comp] = (1.0 - a)*(1.0 - b)*p00[comp] + a*(1.0 - b)*p10[comp]
                          + (1.0 - a)*b*p01[comp] + a*b*p11[comp];
            }
        }
    }

    // --- ZCORN: trilinear resampling per parent cell ----------------------
    const auto parentZ = [&](int i_, int j_, int k_) {
        return zcorn[static_cast<std::size_t>(i_)
                     + 2*static_cast<std::size_t>(nx)*j_
                     + 4*static_cast<std::size_t>(nx)*ny*k_];
    };
    const auto refinedZIndex = [&](int i_, int j_, int k_) {
        return static_cast<std::size_t>(i_)
               + 2*static_cast<std::size_t>(out.dims[0])*j_
               + 4*static_cast<std::size_t>(out.dims[0])*out.dims[1]*k_;
    };

    out.zcorn.resize(8 * static_cast<std::size_t>(out.dims[0]) * out.dims[1] * out.dims[2]);
    for (int kr = 0; kr < out.dims[2]; ++kr) {
        const int ck = request.startIJK[2] + kr / factors[2];
        const int kk = kr % factors[2];
        for (int jr = 0; jr < out.dims[1]; ++jr) {
            const int cj = request.startIJK[1] + jr / factors[1];
            const int jj = jr % factors[1];
            for (int ir = 0; ir < out.dims[0]; ++ir) {
                const int ci = request.startIJK[0] + ir / factors[0];
                const int ii = ir % factors[0];

                for (int dk = 0; dk < 2; ++dk) {
                    const double c = static_cast<double>(kk + dk) / factors[2];
                    for (int dj = 0; dj < 2; ++dj) {
                        const double b = static_cast<double>(jj + dj) / factors[1];
                        for (int di = 0; di < 2; ++di) {
                            const double a = static_cast<double>(ii + di) / factors[0];

                            double z = 0.0;
                            for (int pk = 0; pk < 2; ++pk) {
                                const double wk = (pk == 0) ? (1.0 - c) : c;
                                for (int pj = 0; pj < 2; ++pj) {
                                    const double wj = (pj == 0) ? (1.0 - b) : b;
                                    for (int pi = 0; pi < 2; ++pi) {
                                        const double wi = (pi == 0) ? (1.0 - a) : a;
                                        z += wi*wj*wk*parentZ(2*ci + pi, 2*cj + pj, 2*ck + pk);
                                    }
                                }
                            }
                            out.zcorn[refinedZIndex(2*ir + di, 2*jr + dj, 2*kr + dk)] = z;
                        }
                    }
                }
            }
        }
    }

    // --- ACTNUM: children inherit the parent flag -------------------------
    out.actnum.assign(static_cast<std::size_t>(out.dims[0]) * out.dims[1] * out.dims[2], 1);
    if (actnum) {
        for (int kr = 0; kr < out.dims[2]; ++kr) {
            const int ck = request.startIJK[2] + kr / factors[2];
            for (int jr = 0; jr < out.dims[1]; ++jr) {
                const int cj = request.startIJK[1] + jr / factors[1];
                for (int ir = 0; ir < out.dims[0]; ++ir) {
                    const int ci = request.startIJK[0] + ir / factors[0];
                    const std::size_t parentIdx = static_cast<std::size_t>(ci)
                        + static_cast<std::size_t>(nx)*cj
                        + static_cast<std::size_t>(nx)*ny*ck;
                    const std::size_t childIdx = static_cast<std::size_t>(ir)
                        + static_cast<std::size_t>(out.dims[0])*jr
                        + static_cast<std::size_t>(out.dims[0])*out.dims[1]*kr;
                    out.actnum[childIdx] = actnum[parentIdx];
                }
            }
        }
    }

    return out;
}

} // namespace Refinement
} // namespace Opm
