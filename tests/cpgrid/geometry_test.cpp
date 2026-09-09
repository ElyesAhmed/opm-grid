/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2022 Equinor ASA.

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
#include "config.h"

#define BOOST_TEST_MODULE GeometryTests
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/DefaultGeometryPolicy.hpp>
#include <opm/grid/cpgrid/EntityRep.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>
#include <opm/grid/utility/OpmLog.hpp>

#include <algorithm>

struct Fixture
{
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }

    static int rank()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        return Dune::MPIHelper::instance(m_argc, m_argv).rank();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);
using namespace Dune;

class Null;

BOOST_AUTO_TEST_CASE(vertexgeom)
{
    typedef cpgrid::Geometry<0, 3> Geometry;
    // Default construction.
    [[maybe_unused]] Geometry g_default; // [[maybe_unused]] silences unused warning

    // Construction from point.
    Geometry::GlobalCoordinate c(3.0);
    Geometry g(c);

    // Verification of properties.
    BOOST_CHECK(g.type().isVertex());
    BOOST_CHECK(g.affine());
    BOOST_CHECK_EQUAL(g.corners(), 1);
    BOOST_CHECK_EQUAL(g.corner(0), c);
    Geometry::LocalCoordinate lc(0.0);
    BOOST_CHECK_EQUAL(g.global(lc), c);
    // BOOST_CHECK_THROW(g.local(c), std::exception);
    BOOST_CHECK_EQUAL(g.integrationElement(lc), 1.0);
    BOOST_CHECK_EQUAL(g.volume(), 1.0);
    BOOST_CHECK_EQUAL(g.center(), c);
    // BOOST_CHECK_THROW(g.jacobianTransposed(lc), std::exception);
    // BOOST_CHECK_THROW(g.jacobianInverseTransposed(lc), std::exception);
}


BOOST_AUTO_TEST_CASE(intersectiongeom)
{
    typedef cpgrid::Geometry<2, 3> Geometry;
    // Default construction.
    [[maybe_unused]] Geometry g_default; // [[maybe_unused]] silences unused warning

    // Construction from point and volume.
    Geometry::GlobalCoordinate c(3.0);
    Geometry::ctype v = 8.0;
    Geometry g(c, v);

    // Verification of properties.
    BOOST_CHECK(g.type().isNone());
    BOOST_CHECK(g.affine());
    BOOST_CHECK_EQUAL(g.corners(), 0); // "corners()" returns '0' (None Type -> singular geometry -> no corners)
    // BOOST_CHECK_THROW(g.corner(0), std::exception);
    Geometry::LocalCoordinate lc(0.0);
    BOOST_CHECK_THROW(g.global(lc), std::exception);
    BOOST_CHECK_THROW(g.local(c), std::exception);
    BOOST_CHECK_EQUAL(g.integrationElement(lc), v);
    BOOST_CHECK_EQUAL(g.volume(), v);
    BOOST_CHECK_EQUAL(g.center(), c);
    BOOST_CHECK_THROW(g.jacobianTransposed(lc), std::exception);
    BOOST_CHECK_THROW(g.jacobianInverseTransposed(lc), std::exception);

}


BOOST_AUTO_TEST_CASE(cellgeom)
{
    typedef cpgrid::Geometry<3, 3> Geometry;
    // Default construction.
    Geometry g_default;

    // Construction from point, volume, points and pointindices.
    // Construction from
    // 1. center
    // 2. volume
    // 3. all corners (a pointer of type 'const cpgrid::Geometry<0, 3>')
    // 4. corner indices (a pointer to the first element of "global_refined_cell8corners_indices_storage")
    // First a unit cube, i.e. the mapping represented is the identity.
    typedef Geometry::GlobalCoordinate GC;
    GC c(0.5);
    Geometry::ctype v = 1.0;
    GC corners[8];
    GC cor;
    for (int k = 0; k < 2; ++k) {
        cor[2] = k;
        for (int j = 0; j < 2; ++j) {
            cor[1] = j;
            for (int i = 0; i < 2; ++i) {
                cor[0] = i;
                corners[4*k + 2*j + i] = cor;
            }
        }
    }
//     for (int i = 0; i < 8; ++i) {
//         std::cout << corners[i] << std::endl;
//     }
    // call empty constructor initialze an object that we point to.
    auto pg = std::make_shared<cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3>>();
    (*pg).reserve(8);
    for (const auto& crn : corners)
    {
        (*pg).push_back(cpgrid::Geometry<0, 3>(crn));
    }

    int cor_idx[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    Geometry g(c, v, pg, cor_idx);

    // Verification of properties.
    BOOST_CHECK(g.type().isCube());
    BOOST_CHECK(!g.affine());
    BOOST_CHECK_EQUAL(g.corners(), 8);
    for (int i = 0; i < 8; ++i) {
        BOOST_CHECK_EQUAL(g.corner(i), corners[i]);   // "corner(i)" returns 'pg[cor_idx_[i]].center()'
    }
    BOOST_CHECK_EQUAL(g.volume(), v);
    BOOST_CHECK_EQUAL(g.center(), c);

    // Verification of properties that depend on the mapping.'
    // Construction refined (local) corners.
    typedef Geometry::LocalCoordinate LC;
    const int N = 5;  //  4 in each direction
    const int num_pts = N*N*N;  // total amount of corners
    LC testpts[num_pts] = { LC(0.0) };
    LC pt(0.0);
    for (int k = 0; k < N; ++k) {
        pt[2] = double(k)/double(N-1);
        for (int j = 0; j < N; ++j) {
            pt[1] = double(j)/double(N-1);
            for (int i = 0; i < N; ++i) {
                pt[0] = double(i)/double(N-1);
                testpts[k*N*N + j*N + i] = pt;
//                 std::cout << pt << std::endl;
            }
        }
    }
    Geometry::JacobianTransposed id(0.0);
    id[0][0] = id[1][1] = id[2][2] = 1.0;
    for (int i = 0; i < num_pts; ++i) { // all refined corners
        BOOST_CHECK_EQUAL(g.global(testpts[i]), testpts[i]);
        BOOST_CHECK_EQUAL(g.local(g.global(testpts[i])), testpts[i]);
        BOOST_CHECK_EQUAL(g.integrationElement(testpts[i]), 1.0);   // 'MatrixHelperType::template sqrtDetAAT<3, 3>(testpts[i])'
        BOOST_CHECK_EQUAL(g.jacobianTransposed(testpts[i]), id);
        BOOST_CHECK_EQUAL(g.jacobianInverseTransposed(testpts[i]), id);
    }

    // Next testcase: a degenerate hexahedron, wedge shaped.
    typedef Geometry::GlobalCoordinate GC;
    c = GC(1.0/3.0); c[2] = 0.5;
    v = 0.5;
    corners[5][2] = 0.0;
    corners[7][2] = 0.0;

    // call empty constructor initialze an object that we point to.
    auto pg1 = std::make_shared<cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3>>();
    (*pg1).reserve(8);
    for (const auto& crn : corners)
    {
        (*pg1).push_back(cpgrid::Geometry<0, 3>(crn));
    }
    g = Geometry(c, v, pg1, cor_idx);

    // Verification of properties.
    BOOST_CHECK(g.type().isCube());
    BOOST_CHECK(!g.affine());
    BOOST_CHECK_EQUAL(g.corners(), 8);
    for (int i = 0; i < 8; ++i) {
        BOOST_CHECK_EQUAL(g.corner(i), corners[i]);
    }
    BOOST_CHECK_EQUAL(g.volume(), v);
    BOOST_CHECK_EQUAL(g.center(), c);

    struct Wedge
    {
        static GC global(const LC& lc)
        {
            GC gc(0.0);
            gc[0] = lc[0];
            gc[1] = lc[1];
            gc[2] = (1.0 - lc[0])*lc[2];
            return gc;
        }
        static double integrationElement(const LC& lc)
        {
            return 1.0 - lc[0];
        }
        static Geometry::JacobianTransposed jacobianTransposed(const LC& lc)
        {
            Geometry::JacobianTransposed Jt(0.0);
            Jt[0][0] = 1.0;
            Jt[0][2] = -lc[2];
            Jt[1][1] = 1.0;
            Jt[2][2] = 1.0 - lc[0];
            return Jt;
        }
    };

    // Verification of properties that depend on the mapping.
    const double tolerance = 1e-14;
    for (int i = 0; i < num_pts; ++i) {
        GC gl = g.global(testpts[i]);
        BOOST_CHECK_EQUAL(gl, Wedge::global(testpts[i]));
        BOOST_CHECK_EQUAL(g.integrationElement(testpts[i]), Wedge::integrationElement(testpts[i]));
        Geometry::JacobianTransposed Jt = Wedge::jacobianTransposed(testpts[i]);
        BOOST_CHECK_EQUAL(g.jacobianTransposed(testpts[i]), Jt);
        if (testpts[i][0] < 1.0) {
            // Only do this test if we are away from the degeneracy.
            LC diff = g.local(gl);
            diff -= testpts[i];
            BOOST_CHECK_SMALL(diff.two_norm(), tolerance);
            Geometry::Jacobian Jit = Jt; // This implicitly assumes that the Jacobian is square.
            Jit.invert();
            BOOST_CHECK_EQUAL(g.jacobianInverseTransposed(testpts[i]), Jit);
        }
    }
}
