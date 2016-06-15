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

#if HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif // HAVE_DYNAMIC_BOOST_TEST

#define NVERBOSE

#define BOOST_TEST_MODULE TEST_FLOWDIAGNOSTICSTOOL

#include <boost/test/unit_test.hpp>

#include <opm/flowdiagnostics/FlowDiagnosticsTool.hpp>

#include <opm/flowdiagnostics/CellSet.hpp>
#include <opm/flowdiagnostics/ConnectionValues.hpp>
#include <opm/flowdiagnostics/ConnectivityGraph.hpp>
#include <opm/flowdiagnostics/utility/RandomVector.hpp>

#include <algorithm>

namespace
{
    std::size_t
    numIntConn(const std::size_t nx,
               const std::size_t ny)
    {
        return (nx - 1)*ny + nx*(ny - 1);
    }

    std::vector<int>
    internalConnections(const std::size_t nx,
                        const std::size_t ny)
    {
        auto cellID = [](const std::size_t start,
                         const std::size_t off)
        {
            return static_cast<int>(start + off);
        };

        auto neighbours = std::vector<int>{};
        neighbours.reserve(2 * numIntConn(nx, ny));

        // I connections
        {
            for (auto j = 0*ny; j < ny; ++j) {
                const auto start = j * nx;

                for (auto i = 0*nx + 1; i < nx; ++i) {
                    neighbours.push_back(cellID(start, i - 1));
                    neighbours.push_back(cellID(start, i - 0));
                }
            }
        }

        // J connections
        {
            for (auto j = 0*ny + 1; j < ny; ++j) {
                const auto start = (j - 1)*nx;

                for (auto i = 0*nx; i < nx; ++i) {
                    neighbours.push_back(cellID(start, i + 0 ));
                    neighbours.push_back(cellID(start, i + nx));
                }
            }
        }

        return neighbours;
    }

    std::vector<double>
    flowField(const std::vector<double>::size_type n)
    {
        static Opm::Utility::RandomVector genRandom{};

        return genRandom.normal(n);
    }

} // Namespace anonymous

class Setup
{
public:
    Setup(const std::size_t nx,
          const std::size_t ny);

    const Opm::ConnectivityGraph& connectivity() const;
    const std::vector<double>&    poreVolume()   const;
    const Opm::ConnectionValues&  flux()         const;

private:
    Opm::ConnectivityGraph g_;
    std::vector<double>    pvol_;
    Opm::ConnectionValues  flux_;
};

Setup::Setup(const std::size_t nx,
             const std::size_t ny)
    : g_   (nx * ny, internalConnections(nx, ny))
    , pvol_(g_.numCells(), 0.3)
    , flux_(Opm::ConnectionValues::NumConnections{ g_.numConnections() },
            Opm::ConnectionValues::NumPhases     { 1 })
{
    const auto flux = flowField(g_.numConnections());

    using ConnID = Opm::ConnectionValues::ConnID;

    const auto phaseID =
        Opm::ConnectionValues::PhaseID{ 0 };

    for (decltype(flux_.numConnections())
             conn = 0, nconn = flux_.numConnections();
         conn < nconn; ++conn)
    {
        flux_(ConnID{conn}, phaseID) = flux[conn];
    }
}

const Opm::ConnectivityGraph&
Setup::connectivity() const
{
    return g_;
}

const std::vector<double>&
Setup::poreVolume() const
{
    return pvol_;
}

const Opm::ConnectionValues&
Setup::flux() const
{
    return flux_;
}

BOOST_AUTO_TEST_SUITE(FlowDiagnostics_Tool)

BOOST_AUTO_TEST_CASE (Constructor)
{
    using FDT = Opm::FlowDiagnosticsTool;

    const auto cas = Setup(2, 2);

    FDT diagTool(cas.connectivity());

    diagTool
        .assign(FDT::PoreVolume{ cas.poreVolume() })
        .assign(FDT::ConnectionFlux{ cas.flux() });
}

BOOST_AUTO_TEST_CASE (InjectionDiagnostics)
{
    using FDT = Opm::FlowDiagnosticsTool;

    const auto cas = Setup(2, 2);

    FDT diagTool(cas.connectivity());

    diagTool
        .assign(FDT::PoreVolume{ cas.poreVolume() })
        .assign(FDT::ConnectionFlux{ cas.flux() });

    auto start = std::vector<Opm::CellSet>{};
    {
        start.emplace_back();

        auto& s = start.back();

        s.identify(Opm::CellSetID("I-1"));
        s.insert(0);
    }

    {
        start.emplace_back();

        auto& s = start.back();

        s.identify(Opm::CellSetID("I-2"));
        s.insert(cas.connectivity().numCells() - 1);
    }

    const auto fwd = diagTool
        .computeInjectionDiagnostics(FDT::StartCells{start});

    // Global ToF field (accumulated from all injectors)
    {
        const auto tof = fwd.fd.timeOfFlight();

        BOOST_CHECK_EQUAL(tof.size(), cas.connectivity().numCells());
    }

    // Verify set of start points.
    {
        const auto startpts = fwd.fd.startPoints();

        BOOST_CHECK_EQUAL(startpts.size(), start.size());

        for (const auto& pt : startpts) {
            auto pos =
                std::find_if(start.begin(), start.end(),
                    [&pt](const Opm::CellSet& s)
                    {
                        return s.id().to_string() == pt.to_string();
                    });

            // ID of 'pt' *MUST* be in set of identified start points.
            BOOST_CHECK(pos != start.end());
        }
    }

    // Tracer-ToF
    {
        const auto tof = fwd.fd
            .timeOfFlight(Opm::CellSetID("I-1"));

        for (decltype(tof.cellValueCount())
                 i = 0, n = tof.cellValueCount();
             i < n; ++i)
        {
            const auto v = tof.cellValue(i);

            BOOST_TEST_MESSAGE("[" << i << "] -> ToF["
                               << v.first << "] = "
                               << v.second);
        }
    }

    // Tracer Concentration
    {
        const auto conc = fwd.fd
            .concentration(Opm::CellSetID("I-2"));

        BOOST_TEST_MESSAGE("conc.cellValueCount() = " <<
                           conc.cellValueCount());

        for (decltype(conc.cellValueCount())
                 i = 0, n = conc.cellValueCount();
             i < n; ++i)
        {
            const auto v = conc.cellValue(i);

            BOOST_TEST_MESSAGE("[" << i << "] -> Conc["
                               << v.first << "] = "
                               << v.second);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
