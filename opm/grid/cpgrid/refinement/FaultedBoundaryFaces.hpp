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
#ifndef OPM_GRID_REFINEMENT_FAULTEDBOUNDARYFACES_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_FAULTEDBOUNDARYFACES_HEADER_INCLUDED

#include <array>
#include <vector>

namespace Opm
{
namespace Refinement
{

/// One connection between a refined box-boundary cell and a coarse neighbour
/// cell across a (possibly faulted) box boundary.
struct BoundaryConnection
{
    std::array<int,3> boxCell{};     ///< lattice (i,j,k) in the box's refined frame
    int coarseNeighborCart{-1};      ///< parent Cartesian index of the coarse neighbour
    std::vector<std::array<double,3>> faceNodes; ///< face corner coordinates
};

/// Compute the connections at one box-boundary side using the corner-point
/// processor, so faults on that boundary are handled the same way faults
/// elsewhere are: a refined boundary cell may connect to more than one coarse
/// neighbour (the throw splits the face).
///
/// Builds a mini-grdecl covering the box plus a one-cell coarse shell on the
/// (axis, side) side, refined to cellsPerDim, runs process_grdecl (which splits
/// faces at faults and resolution jumps), and returns the faces between the box
/// boundary cells and the shell (coarse-neighbour) cells. Returns an empty
/// vector if the shell falls outside the parent grid (a domain boundary).
///
/// axis in {0,1,2}; side in {-1, +1}. edgeConformal selects edge-conformal
/// processing, matching the leaf's final edgeConformalizeLeaf pass.
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
                           bool edgeConformal);

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_REFINEMENT_FAULTEDBOUNDARYFACES_HEADER_INCLUDED
