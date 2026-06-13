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

  Tests CpGrid::setPartitionCellGroups - the API that asks load balancing
  to keep groups of cells (e.g. LGR refinement boxes) on one rank, by
  contracting each group into a single partition-graph vertex (PLAN Track
  1 step 6).

  NOTE: full parallel verification (a contracted box staying on one rank
  after scatter) is pending - on synthetic createCartesian grids the
  contraction of a fully-interior region exposes CpGrid's overlap-layer-1
  scatter limitation (see docs/PLAN.md). This test covers the API and the
  no-op (empty groups) path, which must not perturb normal load balancing.
*/
#include <config.h>

#define BOOST_TEST_MODULE PartitionCellGroupsTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <set>
#include <vector>

struct MPIFixture
{
    MPIFixture()
    {
        auto& argv = boost::unit_test::framework::master_test_suite().argv;
        auto& argc = boost::unit_test::framework::master_test_suite().argc;
        Dune::MPIHelper::instance(argc, argv);
    }
};
BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE(roundTrip)
{
    Dune::CpGrid grid;
    grid.createCartesian({{8, 8, 4}}, {{8.0, 8.0, 4.0}});

    BOOST_CHECK(grid.partitionCellGroups().empty());

    std::set<int> box;
    for (int k = 1; k < 3; ++k) {
        for (int j = 2; j < 5; ++j) {
            for (int i = 2; i < 5; ++i) {
                box.insert(i + 8*j + 64*k);
            }
        }
    }
    grid.setPartitionCellGroups({box});
    BOOST_REQUIRE_EQUAL(grid.partitionCellGroups().size(), 1u);
    BOOST_CHECK(grid.partitionCellGroups()[0] == box);
}

BOOST_AUTO_TEST_CASE(emptyGroupsDoNotPerturbLoadBalance)
{
    // With no groups set, load balancing must behave exactly as before:
    // the contraction hook is a no-op. (Regression guard for the
    // GraphOfGridWrappers change.)
    Dune::CpGrid grid;
    grid.createCartesian({{8, 8, 4}}, {{8.0, 8.0, 4.0}});
    BOOST_CHECK(grid.partitionCellGroups().empty());
#if HAVE_MPI
    grid.loadBalance(/*overlapLayers=*/1, Dune::PartitionMethod::zoltanGoG);
#endif
    int local = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        if (element.partitionType() == Dune::InteriorEntity) {
            ++local;
        }
    }
    const int total = grid.comm().sum(local);
    BOOST_CHECK_EQUAL(total, 8*8*4);
}
