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
#include <config.h>

#define BOOST_TEST_MODULE ConformingBuilderTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace
{

struct MPIFixture
{
    MPIFixture()
    {
        auto& argv = boost::unit_test::framework::master_test_suite().argv;
        auto& argc = boost::unit_test::framework::master_test_suite().argc;
        Dune::MPIHelper::instance(argc, argv);
    }
};

struct TestGrdecl
{
    std::array<int,3> dims;
    std::vector<double> coord;
    std::vector<double> zcorn;
    std::vector<int> actnum;

    grdecl raw() const
    {
        grdecl g;
        g.dims[0] = dims[0]; g.dims[1] = dims[1]; g.dims[2] = dims[2];
        g.coord = coord.data();
        g.zcorn = zcorn.data();
        g.actnum = actnum.empty() ? nullptr : actnum.data();
        return g;
    }
};

int cellOf(int doubled) { return doubled / 2; }
int sideOf(int doubled) { return doubled % 2; }

TestGrdecl makeVerticalPillarGrid(const std::array<int,3>& dims,
                                  const std::function<double(int,int,int)>& depth)
{
    TestGrdecl g;
    g.dims = dims;
    const auto& [nx, ny, nz] = dims;

    g.coord.resize(6*(nx + 1)*(ny + 1));
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            double* p = &g.coord[6*(static_cast<std::size_t>(j)*(nx + 1) + i)];
            p[0] = i; p[1] = j; p[2] = 0.0;
            p[3] = i; p[4] = j; p[5] = 100.0;
        }
    }
    g.zcorn.resize(8*static_cast<std::size_t>(nx)*ny*nz);
    for (int k_ = 0; k_ < 2*nz; ++k_) {
        for (int j_ = 0; j_ < 2*ny; ++j_) {
            for (int i_ = 0; i_ < 2*nx; ++i_) {
                g.zcorn[static_cast<std::size_t>(i_) + 2*static_cast<std::size_t>(nx)*j_
                        + 4*static_cast<std::size_t>(nx)*ny*k_] = depth(i_, j_, k_);
            }
        }
    }
    return g;
}

class BuilderGuard
{
public:
    explicit BuilderGuard(std::unique_ptr<Opm::Refinement::Builder> b)
        : previous_{Opm::Refinement::setBuilder(std::move(b))}
    {}
    ~BuilderGuard()
    {
        Opm::Refinement::setBuilder(std::move(previous_));
    }
private:
    std::unique_ptr<Opm::Refinement::Builder> previous_;
};

double totalVolume(const Dune::CpGrid& grid)
{
    double vol = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        vol += element.geometry().volume();
    }
    return vol;
}

// Structural validity of a leaf with a faulted box boundary (phase 1): one
// refined level, volume conserved, every cell closed (sum of area*outward-normal
// is zero), and the box connected to coarse neighbours.
void checkValidFaultedLeaf(const Dune::CpGrid& grid, double volumeBefore)
{
    BOOST_CHECK_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    int refinedToCoarse = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        Dune::FieldVector<double,3> closure(0.0);
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            auto n = is.centerUnitOuterNormal();
            n *= is.geometry().volume();
            closure += n;
            if (is.neighbor() && is.inside().hasFather() != is.outside().hasFather()) {
                ++refinedToCoarse;
            }
        }
        BOOST_CHECK_SMALL(closure.two_norm(), 1e-9);
    }
    BOOST_CHECK_GT(refinedToCoarse, 0);
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE(endToEndSingleBox)
{
    // Distorted vertical-pillar 4x3x3 grid.
    auto depth = [](int i_, int j_, int k_) {
        const double x = cellOf(i_) + sideOf(i_);
        const double y = cellOf(j_) + sideOf(j_);
        return 2.0*(cellOf(k_) + sideOf(k_)) + 0.3*std::sin(0.9*x) + 0.2*std::cos(0.6*y);
    };
    auto parent = makeVerticalPillarGrid({4, 3, 3}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);
    const int cellsBefore = grid.size(0);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // Box of 2x1x1 parents, refined 2x2x2.
    grid.addLgrsUpdateLeafView({{2,2,2}}, {{1,1,1}}, {{3,2,2}}, {"LGR1"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_EQUAL(grid.size(0), cellsBefore - 2 + 2*8);
    BOOST_CHECK_EQUAL(grid.size(0, 0), cellsBefore);    // level 0 view
    BOOST_CHECK_EQUAL(grid.size(1, 0), 16);             // refined level view
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    BOOST_CHECK_EQUAL(grid.getLgrNameToLevel().at("LGR1"), 1);

    // Parent-child structure through the leaf view.
    std::map<int, double> volumePerParent;
    int refinedCount = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        if (element.hasFather()) {
            ++refinedCount;
            const auto father = element.father();
            BOOST_CHECK_EQUAL(father.level(), 0);
            volumePerParent[father.index()] += element.geometry().volume();
            BOOST_CHECK_CLOSE(element.geometryInFather().volume(), 1.0/8.0, 1e-9);
            // Leaf global_cell_ points at the parent's Cartesian index.
            BOOST_CHECK_EQUAL(grid.globalCell()[element.index()],
                              grid.currentData()[0]->globalCell()[father.index()]);
        }
    }
    BOOST_CHECK_EQUAL(refinedCount, 16);
    BOOST_CHECK_EQUAL(volumePerParent.size(), 2u);
    for (const auto& [parentIdx, childVolume] : volumePerParent) {
        const Dune::cpgrid::Entity<0> parentCell(*grid.currentData()[0], parentIdx, true);
        BOOST_CHECK_CLOSE(childVolume, parentCell.geometry().volume(), 1e-8);
    }

    // Intersection consistency over the whole leaf: every interior
    // intersection must be seen from both sides with matching geometry.
    std::map<std::pair<int,int>, int> pairCount;
    double interiorArea = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
            if (intersection.neighbor()) {
                const int inside = intersection.inside().index();
                const int outside = intersection.outside().index();
                BOOST_REQUIRE(inside != outside);
                pairCount[{std::min(inside, outside), std::max(inside, outside)}] += 1;
                interiorArea += intersection.geometry().volume();
            }
        }
    }
    for (const auto& [cells, count] : pairCount) {
        BOOST_CHECK_EQUAL(count, 2); // seen once from each side
    }

    // The coarse i-neighbor of the box (cell (0,1,1)) sees the refined
    // mosaic: 5 unrefined faces + 4 refined ones.
    const auto& dims = grid.logicalCartesianSize();
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        if (element.hasFather()) {
            continue;
        }
        const int cart = grid.globalCell()[element.index()];
        const int i = cart % dims[0];
        const int j = (cart / dims[0]) % dims[1];
        const int k = cart / (dims[0]*dims[1]);
        if (i == 0 && j == 1 && k == 1) {
            int intersections = 0;
            int refinedNeighbors = 0;
            for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
                ++intersections;
                if (intersection.neighbor() && intersection.outside().hasFather()) {
                    ++refinedNeighbors;
                }
            }
            BOOST_CHECK_EQUAL(intersections, 9);
            BOOST_CHECK_EQUAL(refinedNeighbors, 4);
        }
    }

    // Global ids: unique over leaf cells and over leaf points.
    const auto& ids = grid.globalIdSet();
    std::set<std::int64_t> cellIds;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        cellIds.insert(ids.id(element));
    }
    BOOST_CHECK_EQUAL(cellIds.size(), static_cast<std::size_t>(grid.size(0)));

    std::set<std::int64_t> pointIds;
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        pointIds.insert(ids.id(vertex));
    }
    BOOST_CHECK_EQUAL(pointIds.size(), static_cast<std::size_t>(grid.size(3)));
}

BOOST_AUTO_TEST_CASE(stableCellIdDisambiguatesRefinedSiblings)
{
    // global_cell_ is NOT unique on the leaf -- refined siblings share the
    // parent's Cartesian index -- so it cannot key per-cell output gathered
    // across ranks. stableCellId() must give every leaf cell a unique id:
    // coarse cells keep their Cartesian index; refined cells encode
    // (parent Cartesian, child index) under a tag bit.
    auto parent = makeVerticalPillarGrid({4, 3, 3}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));
    grid.addLgrsUpdateLeafView({{2,2,2}}, {{1,1,1}}, {{3,2,2}}, {"LGR1"}); // 2 parents -> 16 refined

    constexpr std::int64_t tag = std::int64_t(1) << 62;
    constexpr std::int64_t childMask = (std::int64_t(1) << 20) - 1;

    const auto sid = grid.stableCellId();
    BOOST_REQUIRE_EQUAL(sid.size(), static_cast<std::size_t>(grid.size(0)));

    // global_cell_ really does collide for siblings, but stableCellId does not.
    const auto& gc = grid.globalCell();
    BOOST_CHECK_LT(std::set<int>(gc.begin(), gc.end()).size(), gc.size());
    BOOST_CHECK_EQUAL(std::set<std::int64_t>(sid.begin(), sid.end()).size(), sid.size());

    std::map<int, std::set<int>> childrenOfParent; // father Cartesian -> child indices seen
    int refined = 0, coarse = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        const int c = element.index();
        if (element.hasFather()) {
            ++refined;
            const int fatherCart = grid.currentData()[0]->globalCell()[element.father().index()];
            BOOST_CHECK((sid[c] & tag) != 0);                       // tagged refined
            BOOST_CHECK_EQUAL((sid[c] & ~tag) >> 20, fatherCart);    // decodes to parent Cartesian
            childrenOfParent[fatherCart].insert(static_cast<int>(sid[c] & childMask));
        }
        else {
            ++coarse;
            BOOST_CHECK_EQUAL(sid[c] & tag, 0);                      // untagged coarse
            BOOST_CHECK_EQUAL(sid[c], static_cast<std::int64_t>(gc[c]));
        }
    }
    BOOST_CHECK_EQUAL(refined, 16);
    BOOST_CHECK_GT(coarse, 0);
    // Each of the two refined parents has all 8 distinct child indices 0..7.
    BOOST_CHECK_EQUAL(childrenOfParent.size(), 2u);
    for (const auto& [father, kids] : childrenOfParent) {
        BOOST_CHECK_EQUAL(kids.size(), 8u);
        BOOST_CHECK_EQUAL(*kids.begin(), 0);
        BOOST_CHECK_EQUAL(*kids.rbegin(), 7);
    }

    // Construction-stable: recomputing yields the identical ids.
    BOOST_CHECK(grid.stableCellId() == sid);
}

BOOST_AUTO_TEST_CASE(twoSeparatedBoxesAndGuards)
{
    auto parent = makeVerticalPillarGrid({6, 3, 3}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // Overlapping boxes are rejected (validation).
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,2,2}, {2,2,2}},
                                                 {{0,0,0}, {1,0,0}},
                                                 {{3,2,2}, {4,2,2}},
                                                 {"A", "B"}),
                      std::invalid_argument);
    BOOST_CHECK_EQUAL(grid.maxLevel(), 0); // grid unchanged

    // Two separated boxes work.
    grid.addLgrsUpdateLeafView({{2,2,2}, {3,3,3}},
                               {{0,0,0}, {4,1,1}},
                               {{1,1,1}, {6,3,3}},
                               {"A", "B"});
    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_EQUAL(grid.getLgrNameToLevel().at("A"), 1);
    BOOST_CHECK_EQUAL(grid.getLgrNameToLevel().at("B"), 2);
    // 54 - 1 - 8 parents + 8 + 8*27 children
    BOOST_CHECK_EQUAL(grid.size(0), 54 - 1 - 8 + 8 + 8*27);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);

    std::set<std::int64_t> cellIds;
    const auto& ids = grid.globalIdSet();
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        cellIds.insert(ids.id(element));
    }
    BOOST_CHECK_EQUAL(cellIds.size(), static_cast<std::size_t>(grid.size(0)));
}

// Two diagonally-adjacent boxes share only an edge (the CARFIN deck
// configuration). The refined corners along the shared edge that are
// interior to a parent cell must merge into single leaf vertices.
BOOST_AUTO_TEST_CASE(edgeSharingBoxesMergeCorners)
{
    auto parent = makeVerticalPillarGrid({4, 4, 2}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // Box A = i,j in [0,2); box B = i,j in [2,4); both k in [0,2).
    // They meet along the vertical edge (i=2, j=2). Factors 2x2x2.
    grid.addLgrsUpdateLeafView({{2,2,2}, {2,2,2}},
                               {{0,0,0}, {2,2,0}},
                               {{2,2,2}, {4,4,2}},
                               {"A", "B"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    // 32 - 2*8 parents + 2*64 children
    BOOST_CHECK_EQUAL(grid.size(0), 32 - 2*8 + 2*64);

    // No two distinct leaf vertices may sit at the same coordinate: the
    // shared-edge interior corners must have merged.
    std::set<std::array<double,3>> coords;
    int vertexCount = 0;
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        const auto& c = vertex.geometry().center();
        coords.insert({c[0], c[1], c[2]});
        ++vertexCount;
    }
    BOOST_CHECK_EQUAL(coords.size(), static_cast<std::size_t>(vertexCount));
    BOOST_CHECK_EQUAL(vertexCount, grid.size(3));

    // Two-sided intersection symmetry over the whole leaf.
    std::map<std::pair<int,int>, int> pairCount;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
            if (intersection.neighbor()) {
                const int in = intersection.inside().index();
                const int out = intersection.outside().index();
                pairCount[{std::min(in, out), std::max(in, out)}] += 1;
            }
        }
    }
    for (const auto& [cells, count] : pairCount) {
        BOOST_CHECK_EQUAL(count, 2);
    }

    // Unique global ids for cells and points.
    const auto& ids = grid.globalIdSet();
    std::set<std::int64_t> cellIds, pointIds;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        cellIds.insert(ids.id(element));
    }
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        pointIds.insert(ids.id(vertex));
    }
    BOOST_CHECK_EQUAL(cellIds.size(), static_cast<std::size_t>(grid.size(0)));
    BOOST_CHECK_EQUAL(pointIds.size(), static_cast<std::size_t>(grid.size(3)));
}

// Face-sharing boxes: the boundary faces of the two refined blocks on the
// shared plane pair into interior faces. The out-of-face subdivision may
// differ between the two boxes; the in-face subdivisions must match.
BOOST_AUTO_TEST_CASE(faceSharingBoxes)
{
    auto parent = makeVerticalPillarGrid({4, 2, 2}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // Box A = i in [0,2), box B = i in [2,4); both full in j,k -> share the
    // i=2 face. Different in-plane-normal factor (rx 2 vs 3), matching
    // in-face factors (ry=rz=2).
    grid.addLgrsUpdateLeafView({{2,2,2}, {3,2,2}},
                               {{0,0,0}, {2,0,0}},
                               {{2,2,2}, {4,2,2}},
                               {"A", "B"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    // Box A (8 parents x 2x2x2) + box B (8 parents x 3x2x2); no coarse left.
    BOOST_CHECK_EQUAL(grid.size(0), 8*8 + 8*12);

    // No coincident-distinct vertices (shared-face corners merged).
    std::set<std::array<double,3>> coords;
    int vertexCount = 0;
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        const auto& c = vertex.geometry().center();
        coords.insert({c[0], c[1], c[2]});
        ++vertexCount;
    }
    BOOST_CHECK_EQUAL(coords.size(), static_cast<std::size_t>(vertexCount));

    // Two-sided intersection symmetry, and the A/B blocks are actually
    // connected (interior faces pairing an A-child with a B-child exist).
    std::map<std::pair<int,int>, int> pairCount;
    int abConnections = 0;
    const auto& dims = grid.logicalCartesianSize();
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
            if (!intersection.neighbor()) {
                continue;
            }
            const int in = intersection.inside().index();
            const int out = intersection.outside().index();
            pairCount[{std::min(in, out), std::max(in, out)}] += 1;
            if (intersection.inside().hasFather() && intersection.outside().hasFather()) {
                const int ci = grid.globalCell()[in] % dims[0];
                const int co = grid.globalCell()[out] % dims[0];
                // One side in box A's i-range [0,2), the other in B's [2,4).
                if ((ci < 2) != (co < 2)) {
                    ++abConnections;
                }
            }
        }
    }
    for (const auto& [cells, count] : pairCount) {
        BOOST_CHECK_EQUAL(count, 2);
    }
    // The shared i=2 face carries ry*rz = 2*2 = 4 connections per parent
    // pair, 2x2 parent pairs in (j,k) -> 16 connections, counted twice.
    BOOST_CHECK_EQUAL(abConnections, 16 * 2);

    const auto& ids = grid.globalIdSet();
    std::set<std::int64_t> cellIds, pointIds;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        cellIds.insert(ids.id(element));
    }
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        pointIds.insert(ids.id(vertex));
    }
    BOOST_CHECK_EQUAL(cellIds.size(), static_cast<std::size_t>(grid.size(0)));
    BOOST_CHECK_EQUAL(pointIds.size(), static_cast<std::size_t>(grid.size(3)));
}

// Non-matching in-face subdivisions are rejected.
BOOST_AUTO_TEST_CASE(faceSharingNonMatchingSubdivisionsThrow)
{
    auto parent = makeVerticalPillarGrid({4, 2, 2}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // ry differs (2 vs 3) on the shared i=2 face.
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,2,2}, {2,3,2}},
                                                 {{0,0,0}, {2,0,0}},
                                                 {{2,2,2}, {4,2,2}},
                                                 {"A", "B"}),
                      std::logic_error);
    BOOST_CHECK_EQUAL(grid.maxLevel(), 0);
}

// Nested LGR (a box whose parent is another LGR). The level grids are built
// over their parent (Phase A/B of docs/NESTED_LGR_PLAN.md), but the recursive
// leaf stitching (Phase C) is not implemented yet, so the build is expected to
// reach that boundary and throw a precise, staged message. Locks in: (a) parent
// ordering is honoured, (b) the staged failure is the *leaf* one (i.e. the
// nested level grids were successfully assembled first).
BOOST_AUTO_TEST_CASE(nestedRefinementReachesLeafBoundary)
{
    auto parent = makeVerticalPillarGrid({4, 2, 2}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // LGR1: 2x1x1 GLOBAL parents -> local refined dims 4x2x2.
    // NEST1: a box inside LGR1's local space, parent "LGR1".
    const auto isNestedLeafBoundary = [](const std::logic_error& e) {
        return std::string(e.what()).find("Nested LGR leaf assembly is not implemented")
               != std::string::npos;
    };
    BOOST_CHECK_EXCEPTION(
        grid.addLgrsUpdateLeafView({{2,2,2}, {2,2,2}},
                                   {{1,0,0}, {0,0,0}},
                                   {{3,1,1}, {2,2,2}},
                                   {"LGR1", "NEST1"},
                                   {"GLOBAL", "LGR1"}),
        std::logic_error, isNestedLeafBoundary);
}

// The builder appends level grids in request order and resolves each box's
// parent by name, so a child must follow its parent. A child listed before its
// parent must be rejected with an actionable message (the simulator's
// topological sort is what guarantees this never happens in practice).
BOOST_AUTO_TEST_CASE(nestedChildBeforeParentThrows)
{
    auto parent = makeVerticalPillarGrid({4, 2, 2}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    const auto isUnbuiltParent = [](const std::logic_error& e) {
        return std::string(e.what()).find("has not been built yet")
               != std::string::npos;
    };
    BOOST_CHECK_EXCEPTION(
        grid.addLgrsUpdateLeafView({{2,2,2}, {2,2,2}},
                                   {{0,0,0}, {1,0,0}},
                                   {{2,2,2}, {3,1,1}},
                                   {"NEST1", "LGR1"},
                                   {"LGR1", "GLOBAL"}),
        std::logic_error, isUnbuiltParent);
}

BOOST_AUTO_TEST_CASE(faultInsideBoxEndToEnd)
{
    // Fault between i=1 and i=2 columns, throw 0.6; box covers it.
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto parent = makeVerticalPillarGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // The box {1..3} has the fault *inside* (between i=1 and i=2) and
    // unfaulted boundaries (i=0|1 and i=3|4 columns are unfaulted... note
    // the fault plane is i=2, the box boundary planes are i=1 and i=3,
    // both unfaulted since the throw is constant on each side).
    grid.addLgrsUpdateLeafView({{2,2,2}}, {{1,0,0}}, {{3,2,2}}, {"LGR1"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);

    // Both-sided intersection consistency across the refined fault.
    std::map<std::pair<int,int>, int> pairCount;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
            if (intersection.neighbor()) {
                const int inside = intersection.inside().index();
                const int outside = intersection.outside().index();
                pairCount[{std::min(inside, outside), std::max(inside, outside)}] += 1;
            }
        }
    }
    for (const auto& [cells, count] : pairCount) {
        BOOST_CHECK_EQUAL(count, 2);
    }
}

BOOST_AUTO_TEST_CASE(faultAtBoxBoundaryBuilds)
{
    // Fault between i=1 and i=2 (throw 0.6); box boundary right on the fault
    // plane. The boundary parents each connect to two coarse cells, so the
    // boundary is rebuilt as split faces. Phase 1: assert the leaf is a valid
    // polyhedral grid (builds, volume conserved, every interior face two-sided,
    // refined boundary cells connected to the coarse neighbours).
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto parent = makeVerticalPillarGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{2,0,0}}, {{4,2,2}}, {"LGR1"}));
    BOOST_CHECK_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);

    // Geometric validity: every leaf cell is closed, i.e. the sum of
    // area * outward-unit-normal over all its faces is zero (divergence
    // theorem). This validates the synthetic split faces' areas AND normal
    // orientation, independent of how the interface is subdivided. Also: every
    // interior face is two-sided, and the box connects to coarse neighbours.
    std::map<std::pair<int,int>, int> pairCount;
    int refinedToCoarse = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        Dune::FieldVector<double,3> closure(0.0);
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            auto n = is.centerUnitOuterNormal();
            n *= is.geometry().volume();
            closure += n;
            if (is.neighbor()) {
                const int a = is.inside().index();
                const int b = is.outside().index();
                BOOST_REQUIRE(a != b);
                pairCount[{std::min(a, b), std::max(a, b)}] += 1;
                if (is.inside().hasFather() != is.outside().hasFather()) {
                    ++refinedToCoarse;
                }
            }
        }
        BOOST_CHECK_SMALL(closure.two_norm(), 1e-9);
    }
    for (const auto& [cells, count] : pairCount) {
        BOOST_CHECK_EQUAL(count % 2, 0); // each shared face seen from both sides
    }
    BOOST_CHECK_GT(refinedToCoarse, 0);  // the box connects to coarse neighbours
}

// Note on the boundary throw branches: with vertical pillars + a ZCORN throw
// (the helper here), a faulted box boundary always produces *split* faces (a
// parent connecting to two coarse cells) — the 'fault-split face' branch. The
// sibling 'partial (faulted or degenerate) face' branch needs the geometric
// offset-pillar fault that the helper cannot express; it is exercised by the
// flow deck opm-tests/flow_diagnostic_test/SIMPLE_2PH_W_FAULT_LGR.DATA.

BOOST_AUTO_TEST_CASE(faultAtBoxBoundaryJDirectionBuilds)
{
    // Fault and box boundary in the J direction — the other lateral axis.
    auto depth = [](int i_, int j_, int k_) {
        (void)i_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(j_) >= 2) ? base + 0.6 : base;
    };
    auto parent = makeVerticalPillarGrid({2, 4, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{0,2,0}}, {{2,4,2}}, {"LGR1"}));
    checkValidFaultedLeaf(grid, volumeBefore);
}

BOOST_AUTO_TEST_CASE(faultAtBoxBoundaryUpThrownSideBuilds)
{
    // Box on the *up-thrown* side: fault between i=1 and i=2 (cells i>=2 thrown
    // down), box i0-1 so its right boundary (i=2) faces the thrown neighbour.
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto parent = makeVerticalPillarGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{0,0,0}}, {{2,2,2}}, {"LGR1"}));
    checkValidFaultedLeaf(grid, volumeBefore);
}

BOOST_AUTO_TEST_CASE(faultLargeThrowAtBoxBoundaryBuilds)
{
    // Throw larger than one layer (3.0 > layer thickness 2.0), three layers, so
    // a boundary parent spans two coarse neighbours — a wider overlap pattern.
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 3.0 : base;
    };
    auto parent = makeVerticalPillarGrid({4, 2, 3}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{2,0,0}}, {{4,2,3}}, {"LGR1"}));
    checkValidFaultedLeaf(grid, volumeBefore);
}

BOOST_AUTO_TEST_CASE(faultNotOnBoxBoundaryBuilds)
{
    // Control: a fault exists (between i=0 and i=1) but the box (i=2..3) sits
    // entirely on the down-thrown side, so its boundary toward the coarse
    // neighbour (i=1) is a full, matching face. This must build today and
    // conserve volume — it isolates "fault on the boundary" from "fault merely
    // present in the grid".
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 1) ? base + 0.6 : base;
    };
    auto parent = makeVerticalPillarGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    BOOST_CHECK_NO_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{2,0,0}}, {{4,2,2}}, {"LGR1"}));
    BOOST_CHECK_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
}

BOOST_AUTO_TEST_CASE(processGrdeclSplitsFaultedCoarseFineColumnPair)
{
    // Kernel proof for faults-at-box-boundary (docs/PLAN.md Track 1): the planned
    // fix builds the connecting faces at a faulted box boundary by running
    // process_grdecl on a mini-grdecl (coarse neighbour + refined box boundary)
    // and reading the split connections. This checks that mechanism in isolation.
    //
    // Column 0 (coarse): one cell z[0,4]; its second logical cell is pinched
    //   (z[4,4]) so process_grdecl removes it -> one coarse cell.
    // Column 1 (refined): two cells z[1,3] and z[3,5] -- the same span shifted by
    //   a 1.0 throw. The coarse cell must connect to BOTH refined cells.
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        if (cellOf(i_) == 0) {                 // coarse column: cell0 [0,4], cell1 pinched at 4
            return (k_ == 0) ? 0.0 : 4.0;
        }
        static const double fine[4] = {1.0, 3.0, 3.0, 5.0};   // refined, throw +1
        return fine[k_];
    };
    auto g = makeVerticalPillarGrid({2, 1, 2}, depth);
    auto raw = g.raw();

    struct processed_grid out;
    const int ok = process_grdecl(/*pinchActive=*/1, /*edge_conformal=*/0, 1e-6, &raw, nullptr, &out);
    BOOST_REQUIRE(ok);

    // Pinched coarse cell removed: one coarse + two refined.
    BOOST_CHECK_EQUAL(out.number_of_cells, 3);

    // The split connection lives in the interior I-faces between the two columns
    // (face_neighbors holds compressed active cell ids; -1 is a domain boundary).
    // The coarse cell (0) must connect to both refined cells (1 and 2).
    int interiorIFaces = 0;
    std::set<int> coarseConnectedTo;
    for (unsigned f = 0; f < out.number_of_faces; ++f) {
        const int a = out.face_neighbors[2*f];
        const int b = out.face_neighbors[2*f + 1];
        if (out.face_tag[f] == I_FACE && a >= 0 && b >= 0) {
            ++interiorIFaces;
            BOOST_CHECK((a == 0) || (b == 0));        // one side is the coarse cell
            coarseConnectedTo.insert(a == 0 ? b : a);
        }
    }
    BOOST_CHECK_EQUAL(interiorIFaces, 2);
    BOOST_CHECK((coarseConnectedTo == std::set<int>{1, 2}));

    free_processed_grid(&out);
}
