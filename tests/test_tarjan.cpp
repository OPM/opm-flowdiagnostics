/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define NVERBOSE

#define BOOST_TEST_MODULE TarjanImplementationTest

#include <boost/test/unit_test.hpp>

#include <opm/flowdiagnostics/reorder/tarjan.h>

#include <memory>

namespace {
    struct DestroyWorkSpace
    {
        void operator()(TarjanWorkSpace* ws);
    };

    struct DestroySCCResult
    {
        void operator()(TarjanSCCResult *scc);
    };

    void DestroyWorkSpace::operator()(TarjanWorkSpace* ws)
    {
        destroy_tarjan_workspace(ws);
    }

    void DestroySCCResult::operator()(TarjanSCCResult* scc)
    {
        destroy_tarjan_sccresult(scc);
    }

    using WorkSpace =
        std::unique_ptr<TarjanWorkSpace, DestroyWorkSpace>;

    using SCCResult =
        std::unique_ptr<TarjanSCCResult, DestroySCCResult>;
} // Anonymous

BOOST_AUTO_TEST_SUITE(Two_By_Two)

// +-----+-----+
// |  2  |  3  |
// +-----+-----+
// |  0  |  1  |
// +-----+-----+

BOOST_AUTO_TEST_CASE (FullySeparable)
{
    // Quarter five-spot pattern:
    //   0 -> 1
    //   0 -> 2
    //   1 -> 3
    //   2 -> 3
    //
    const int ia[] = { 0, 0, 1, 2, 4 };
    const int ja[] = { 0, 0, 1, 2 };

    const std::size_t nv           = (sizeof ia) / (sizeof ia[0]) - 1;
    const std::size_t expect_ncomp = 4;

    auto scc = SCCResult{ tarjan(nv, ia, ja) };

    const auto ncomp = tarjan_get_numcomponents(scc.get());

    BOOST_CHECK_EQUAL(ncomp, expect_ncomp);

    const int expect_vert[] = { 0, 1, 2, 3 };

    for (auto comp = 0*ncomp; comp < ncomp; ++comp) {
        const auto c = tarjan_get_strongcomponent(scc.get(), comp);

        BOOST_CHECK_EQUAL(c.size     , 1);
        BOOST_CHECK_EQUAL(c.vertex[0], expect_vert[comp]);
    }
}

BOOST_AUTO_TEST_CASE (Loop)
{
    // Circulation:
    //   0 -> 1
    //   1 -> 3
    //   3 -> 2
    //   2 -> 0
    //
    const int ia[] = { 0, 1, 2, 3, 4 };
    const int ja[] = { 2, 0, 3, 1 };

    const std::size_t nv           = (sizeof ia) / (sizeof ia[0]) - 1;
    const std::size_t expect_ncomp = 1;

    auto scc = SCCResult{ tarjan(nv, ia, ja) };

    const auto ncomp = tarjan_get_numcomponents(scc.get());

    BOOST_CHECK_EQUAL(ncomp, expect_ncomp);

    const auto c = tarjan_get_strongcomponent(scc.get(), 0);

    BOOST_CHECK_EQUAL(c.size, nv);

    // Cell indices within component returned in (essentially) arbitrary
    // order.  This particular order happened to be correct at the time the
    // test was implemented so the assertion on 'vertex' is only usable as a
    // regression test.
    const int expect_vert[] = { 1, 3, 2, 0 };
    BOOST_CHECK_EQUAL_COLLECTIONS(c.vertex   , c.vertex + c.size,
                                  expect_vert, expect_vert + nv);
}

BOOST_AUTO_TEST_CASE (DualPath)
{
    // Two flow paths from cell 0 to cell 2:
    //   0 -> 1
    //   1 -> 3
    //   3 -> 2
    //   0 -> 2
    //
    const int ia[] = { 0, 0, 2, 3, 4 };
    const int ja[] = { 0, 0, 3, 1 };

    const std::size_t nv           = (sizeof ia) / (sizeof ia[0]) - 1;
    const std::size_t expect_ncomp = 4;

    auto scc = SCCResult{ tarjan(nv, ia, ja) };

    const auto ncomp = tarjan_get_numcomponents(scc.get());

    BOOST_CHECK_EQUAL(ncomp, expect_ncomp);

    // Cell 0 is a source and cell 2 is a sink so first and last cells must
    // be 0 and 2 respectively.  The order of cells 1 and 3 is determined by
    // the flow path in which 1 precedes 3.

    const int expect_vert[] = { 0, 1, 3, 2 };

    for (auto comp = 0*ncomp; comp < ncomp; ++comp) {
        const auto c = tarjan_get_strongcomponent(scc.get(), comp);

        BOOST_CHECK_EQUAL(c.size     , 1);
        BOOST_CHECK_EQUAL(c.vertex[0], expect_vert[comp]);
    }
}

BOOST_AUTO_TEST_CASE (IsolatedFlows)
{
    // Compartmentalised reservoir with source and sink in each compartment.
    //   0 -> 2
    //   3 -> 1
    //
    const int ia[] = { 0, 0, 1, 2, 2 };
    const int ja[] = { 3, 0 };

    const std::size_t nv           = (sizeof ia) / (sizeof ia[0]) - 1;
    const std::size_t expect_ncomp = 4;

    auto scc = SCCResult{ tarjan(nv, ia, ja) };

    const auto ncomp = tarjan_get_numcomponents(scc.get());

    BOOST_CHECK_EQUAL(ncomp, expect_ncomp);

    // Sources before sinks, but no a priori ordering between sources or
    // between sinks.  This particular order happened to be correct at the
    // time the test was implemented so the assertion on 'vert' is only
    // usable as a regression test.
    const int expect_vert[] = { 0, 3, 1, 2 };

    for (auto comp = 0*ncomp; comp < ncomp; ++comp) {
        const auto c = tarjan_get_strongcomponent(scc.get(), comp);

        BOOST_CHECK_EQUAL(c.size     , 1);
        BOOST_CHECK_EQUAL(c.vertex[0], expect_vert[comp]);
    }
}

BOOST_AUTO_TEST_SUITE_END()
