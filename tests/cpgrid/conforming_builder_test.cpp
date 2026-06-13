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

    // Touching boxes are rejected by the builder.
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,2,2}, {2,2,2}},
                                                 {{0,0,0}, {2,0,0}},
                                                 {{2,2,2}, {4,2,2}},
                                                 {"A", "B"}),
                      std::logic_error);
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

// Face-sharing boxes are rejected (mosaic-mosaic pairing not implemented).
BOOST_AUTO_TEST_CASE(faceSharingBoxesThrow)
{
    auto parent = makeVerticalPillarGrid({4, 2, 2}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // Box A = i in [0,2), box B = i in [2,4); both full in j,k -> share the
    // i=2 face.
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,2,2}, {2,2,2}},
                                                 {{0,0,0}, {2,0,0}},
                                                 {{2,2,2}, {4,2,2}},
                                                 {"A", "B"}),
                      std::logic_error);
    BOOST_CHECK_EQUAL(grid.maxLevel(), 0);
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

BOOST_AUTO_TEST_CASE(faultAtBoxBoundaryThrows)
{
    // Fault between i=1 and i=2; box boundary right on the fault plane.
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto parent = makeVerticalPillarGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    BuilderGuard guard(std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{2,0,0}}, {{4,2,2}}, {"LGR1"}),
                      std::logic_error);
}
