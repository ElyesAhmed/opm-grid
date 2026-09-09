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

#define BOOST_TEST_MODULE RefinementSeamTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <memory>
#include <stdexcept>
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

// Builder that only records what it was asked to do.
class RecordingBuilder : public Opm::Refinement::Builder
{
public:
    void build(Dune::CpGrid&, const std::vector<Opm::Refinement::BlockRefinement>& requests) override
    {
        recorded = requests;
        ++callCount;
    }

    std::vector<Opm::Refinement::BlockRefinement> recorded{};
    int callCount = 0;
};

// RAII registration so a test cannot leak its builder into the next one.
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

Dune::CpGrid makeGrid()
{
    Dune::CpGrid grid;
    grid.createCartesian({4, 3, 3}, {1.0, 1.0, 1.0});
    return grid;
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE(throwsWithoutRegisteredBuilder)
{
    auto grid = makeGrid();
    BOOST_CHECK(Opm::Refinement::builder() == nullptr);
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{1,1,1}}, {{3,2,2}}, {"LGR1"}),
                      std::logic_error);
    BOOST_CHECK_EQUAL(grid.maxLevel(), 0);
}

BOOST_AUTO_TEST_CASE(forwardsParsedRequestsToBuilder)
{
    auto grid = makeGrid();
    auto recorder = std::make_unique<RecordingBuilder>();
    auto* recorderPtr = recorder.get();
    BuilderGuard guard(std::move(recorder));

    grid.addLgrsUpdateLeafView({{2,2,2}, {3,3,3}},
                               {{0,0,0}, {2,0,0}},
                               {{1,1,1}, {4,2,2}},
                               {"LGR1", "LGR2"});

    BOOST_REQUIRE_EQUAL(recorderPtr->callCount, 1);
    BOOST_REQUIRE_EQUAL(recorderPtr->recorded.size(), 2u);
    BOOST_CHECK_EQUAL(recorderPtr->recorded[0].name, "LGR1");
    BOOST_CHECK_EQUAL(recorderPtr->recorded[0].parentGridName, "GLOBAL");
    BOOST_CHECK_EQUAL(recorderPtr->recorded[1].name, "LGR2");
    BOOST_CHECK(recorderPtr->recorded[1].cellsPerDim == (std::array<int,3>{3,3,3}));
    BOOST_CHECK(recorderPtr->recorded[1].startIJK == (std::array<int,3>{2,0,0}));
    BOOST_CHECK(recorderPtr->recorded[1].endIJK == (std::array<int,3>{4,2,2}));
}

BOOST_AUTO_TEST_CASE(rejectsOverlappingBoxesBeforeCallingBuilder)
{
    auto grid = makeGrid();
    auto recorder = std::make_unique<RecordingBuilder>();
    auto* recorderPtr = recorder.get();
    BuilderGuard guard(std::move(recorder));

    // Boxes [0,2)x[0,2)x[0,2) and [1,3)x[1,2)x[1,2) overlap.
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,2,2}, {2,2,2}},
                                                 {{0,0,0}, {1,1,1}},
                                                 {{2,2,2}, {3,2,2}},
                                                 {"LGR1", "LGR2"}),
                      std::invalid_argument);
    BOOST_CHECK_EQUAL(recorderPtr->callCount, 0);
}

BOOST_AUTO_TEST_CASE(rejectsInvalidBoxesAndSizes)
{
    auto grid = makeGrid();
    auto recorder = std::make_unique<RecordingBuilder>();
    BuilderGuard guard(std::move(recorder));

    // start == end
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{1,1,1}}, {{1,2,2}}, {"LGR1"}),
                      std::invalid_argument);
    // non-positive subdivisions
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{0,2,2}}, {{1,1,1}}, {{2,2,2}}, {"LGR1"}),
                      std::invalid_argument);
    // mismatching vector sizes
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{1,1,1}, {0,0,0}}, {{2,2,2}}, {"LGR1"}),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(disjointBoxesTouchingAtFaceAreAccepted)
{
    auto grid = makeGrid();
    auto recorder = std::make_unique<RecordingBuilder>();
    auto* recorderPtr = recorder.get();
    BuilderGuard guard(std::move(recorder));

    // [0,2) and [2,4) in I touch at the face I=2 but do not overlap.
    grid.addLgrsUpdateLeafView({{2,2,2}, {2,2,2}},
                               {{0,0,0}, {2,0,0}},
                               {{2,2,2}, {4,2,2}},
                               {"LGR1", "LGR2"});
    BOOST_CHECK_EQUAL(recorderPtr->callCount, 1);
}
