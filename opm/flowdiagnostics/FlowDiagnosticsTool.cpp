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

#include <opm/flowdiagnostics/FlowDiagnosticsTool.hpp>

#include <opm/flowdiagnostics/CellSet.hpp>
#include <opm/flowdiagnostics/ConnectionValues.hpp>
#include <opm/flowdiagnostics/ConnectivityGraph.hpp>

#include <algorithm>
#include <exception>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <opm/flowdiagnostics/utility/RandomVector.hpp>

namespace { namespace Mock {

    std::vector<double>
    FieldValue(const std::vector<double>::size_type n,
               const double mean, const double stdev)
    {
        static Opm::Utility::RandomVector genRandom{};

        return genRandom.normal(n, mean, stdev);
    }

    std::vector<int>
    Index(const std::vector<double>::size_type n,
          const int maxIdx)
    {
        static Opm::Utility::RandomVector genRandom{};

        return genRandom.index(n, maxIdx);
    }

}} // namespace (anonymous)::Mock

// ---------------------------------------------------------------------
// Class Opm::FlowDiagnosticsSolution::Impl
// ---------------------------------------------------------------------

class Opm::FlowDiagnosticsSolution::Impl
{
public:
    using GlobalToF = std::vector<double>;

    struct TimeOfFlight
    {
        CellSetValues data;
    };

    struct Concentration
    {
        CellSetValues data;
    };

    void assignToF(GlobalToF&& tof);
    void assign   (const CellSetID& i, TimeOfFlight&&  tof);
    void assign   (const CellSetID& i, Concentration&& conc);

    std::vector<CellSetID> startPoints() const;

    const GlobalToF& timeOfFlight() const;

    CellSetValues timeOfFlight (const CellSetID& tracer) const;
    CellSetValues concentration(const CellSetID& tracer) const;

private:
    struct CompareCellSetIDs
    {
        bool operator()(const CellSetID& x,
                        const CellSetID& y) const
        {
            return x.to_string() < y.to_string();
        }
    };

    using SolutionMap =
        std::map<CellSetID, CellSetValues, CompareCellSetIDs>;

    GlobalToF   tof_;
    SolutionMap tracerToF_;
    SolutionMap tracer_;

    void assign(const CellSetID& i,
                CellSetValues&&  x,
                SolutionMap&     soln);

    CellSetValues
    solutionValues(const CellSetID&   i,
                   const SolutionMap& soln) const;
};

void
Opm::FlowDiagnosticsSolution::Impl::assignToF(GlobalToF&& tof)
{
    tof_ = std::move(tof);
}

void
Opm::FlowDiagnosticsSolution::
Impl::assign(const CellSetID& i, TimeOfFlight&& tof)
{
    assign(i, std::move(tof.data), tracerToF_);
}

void
Opm::FlowDiagnosticsSolution::
Impl::assign(const CellSetID& i, Concentration&& conc)
{
    assign(i, std::move(conc.data), tracer_);
}

std::vector<Opm::CellSetID>
Opm::FlowDiagnosticsSolution::Impl::startPoints() const
{
    auto s = std::vector<CellSetID>{};
    s.reserve(tracer_.size());

    for (const auto& t : tracer_) {
        s.emplace_back(t.first);
    }

    return s;
}

const Opm::FlowDiagnosticsSolution::Impl::GlobalToF&
Opm::FlowDiagnosticsSolution::Impl::timeOfFlight() const
{
    return tof_;
}

Opm::CellSetValues
Opm::FlowDiagnosticsSolution::
Impl::timeOfFlight(const CellSetID& tracer) const
{
    return solutionValues(tracer, tracerToF_);
}

Opm::CellSetValues
Opm::FlowDiagnosticsSolution::
Impl::concentration(const CellSetID& tracer) const
{
    return solutionValues(tracer, tracer_);
}

void
Opm::FlowDiagnosticsSolution::
Impl::assign(const CellSetID& i,
             CellSetValues&&  x,
             SolutionMap&     soln)
{
    soln[i] = std::move(x);
}

Opm::CellSetValues
Opm::FlowDiagnosticsSolution::
Impl::solutionValues(const CellSetID&   i,
                     const SolutionMap& soln) const
{
    auto p = soln.find(i);

    if (p == soln.end()) {
        return CellSetValues{};
    }

    return p->second;
}

// ---------------------------------------------------------------------
// Class Opm::FlowDiagnosticsTool::Impl
// ---------------------------------------------------------------------

class Opm::FlowDiagnosticsTool::Impl
{
public:
    explicit Impl(ConnectivityGraph g);

    void assign(const PoreVolume&     pv);
    void assign(const ConnectionFlux& flux);

    Forward injDiag (const StartCells& start);
    Reverse prodDiag(const StartCells& start);

private:
    ConnectivityGraph g_;

    std::vector<double> pvol_;
    ConnectionValues    flux_;
};

Opm::FlowDiagnosticsTool::Impl::Impl(ConnectivityGraph g)
    : g_   (std::move(g))
    , pvol_(g_.numCells(), 0.0)
    , flux_(ConnectionValues::NumConnections{ 0 },
            ConnectionValues::NumPhases     { 0 })
{}

void
Opm::FlowDiagnosticsTool::Impl::assign(const PoreVolume& pv)
{
    if (pv.data.size() != pvol_.size()) {
        throw std::logic_error("Inconsistently sized input "
                               "pore-volume field");
    }

    pvol_ = pv.data;
}

void
Opm::FlowDiagnosticsTool::Impl::assign(const ConnectionFlux& flux)
{
    if (flux.data.numConnections() != g_.numConnections()) {
        throw std::logic_error("Inconsistently sized input "
                               "flux field");
    }

    flux_ = flux.data;
}

Opm::FlowDiagnosticsTool::Forward
Opm::FlowDiagnosticsTool::Impl::injDiag(const StartCells& start)
{
    using SampleSize = Opm::Utility::RandomVector::Size;
    using Soln       = FlowDiagnosticsSolution::Impl;
    using ToF        = Soln::TimeOfFlight;
    using Conc       = Soln::Concentration;

    using SolnPtr = std::unique_ptr<Soln>;

    SolnPtr x(new Soln());

    const auto avgToF = 20.0e3;
    const auto stdToF = 500.0;

    x->assignToF(Mock::FieldValue(g_.numCells(), avgToF, stdToF));

    for (const auto& pt : start.points)
    {
        const auto npts = static_cast<SampleSize>
            (std::distance(pt.begin(), pt.end()));

        const auto n = std::min(
            { npts, g_.numCells(), static_cast<SampleSize>(100) }
        );

        const auto idx = Mock::Index(n, g_.numCells() - 1);

        {
            const auto tof = Mock::FieldValue(n, avgToF, stdToF);

            auto val = CellSetValues{ n };
            for (auto i = 0*n; i < n; ++i) {
                val.addCellValue(idx[i], tof[i]);
            }

            x->assign(pt.id(), ToF{ val });
        }

        {
            const auto conc = Mock::FieldValue(n, 0.5, 0.15);

            auto val = CellSetValues{ n };
            for (auto i = 0*n; i < n; ++i) {
                val.addCellValue(idx[i], conc[i]);
            }

            x->assign(pt.id(), Conc{ val });
        }
    }

    return Forward{ FlowDiagnosticsSolution(std::move(x)) };
}

Opm::FlowDiagnosticsTool::Reverse
Opm::FlowDiagnosticsTool::Impl::prodDiag(const StartCells& start)
{
    using SampleSize = Opm::Utility::RandomVector::Size;
    using Soln       = FlowDiagnosticsSolution::Impl;
    using ToF        = Soln::TimeOfFlight;
    using Conc       = Soln::Concentration;

    using SolnPtr = std::unique_ptr<Soln>;

    SolnPtr x(new Soln());

    const auto avgToF = 20.0e3;
    const auto stdToF = 500.0;

    x->assignToF(Mock::FieldValue(g_.numCells(), avgToF, stdToF));

    for (const auto& pt : start.points)
    {
        const auto npts = static_cast<SampleSize>
            (std::distance(pt.begin(), pt.end()));

        const auto n = std::min(
            { npts, g_.numCells(), static_cast<SampleSize>(100) }
        );

        const auto idx = Mock::Index(n, g_.numCells() - 1);

        {
            const auto tof = Mock::FieldValue(n, avgToF, stdToF);

            auto val = CellSetValues{ n };
            for (auto i = 0*n; i < n; ++i) {
                val.addCellValue(idx[i], tof[i]);
            }

            x->assign(pt.id(), ToF{ val });
        }

        {
            const auto conc = Mock::FieldValue(n, 0.5, 0.15);

            auto val = CellSetValues{ n };
            for (auto i = 0*n; i < n; ++i) {
                val.addCellValue(idx[i], conc[i]);
            }

            x->assign(pt.id(), Conc{ val });
        }
    }

    return Reverse{ FlowDiagnosticsSolution(std::move(x)) };
}

// =====================================================================
// Implementation of public interface below separator
// =====================================================================

// ---------------------------------------------------------------------
// Class Opm::FlowDiagnosticsSolution
// ---------------------------------------------------------------------

Opm::FlowDiagnosticsSolution::~FlowDiagnosticsSolution()
{}

Opm::FlowDiagnosticsSolution::
FlowDiagnosticsSolution(const FlowDiagnosticsSolution& rhs)
    : pImpl_(new Impl(*rhs.pImpl_))
{}

Opm::FlowDiagnosticsSolution::
FlowDiagnosticsSolution(FlowDiagnosticsSolution&& rhs)
    : pImpl_(std::move(rhs.pImpl_))
{}

Opm::FlowDiagnosticsSolution::
FlowDiagnosticsSolution(std::unique_ptr<Impl> pImpl)
    : pImpl_(std::move(pImpl))
{}

std::vector<Opm::CellSetID>
Opm::FlowDiagnosticsSolution::startPoints() const
{
    return pImpl_->startPoints();
}

const std::vector<double>&
Opm::FlowDiagnosticsSolution::timeOfFlight() const
{
    return pImpl_->timeOfFlight();
}

Opm::CellSetValues
Opm::FlowDiagnosticsSolution::timeOfFlight(const CellSetID& tracer) const
{
    return pImpl_->timeOfFlight(tracer);
}

Opm::CellSetValues
Opm::FlowDiagnosticsSolution::concentration(const CellSetID& tracer) const
{
    return pImpl_->concentration(tracer);
}

// ---------------------------------------------------------------------
// Class Opm::FlowDiagnosticsTool
// ---------------------------------------------------------------------

Opm::FlowDiagnosticsTool::
FlowDiagnosticsTool(const ConnectivityGraph& conn)
    : pImpl_(new Impl(conn))
{}

Opm::FlowDiagnosticsTool::~FlowDiagnosticsTool()
{}

Opm::FlowDiagnosticsTool&
Opm::FlowDiagnosticsTool::assign(const PoreVolume& pv)
{
    pImpl_->assign(pv);

    return *this;
}

Opm::FlowDiagnosticsTool&
Opm::FlowDiagnosticsTool::assign(const ConnectionFlux& flux)
{
    pImpl_->assign(flux);

    return *this;
}

Opm::FlowDiagnosticsTool::Forward
Opm::FlowDiagnosticsTool::
computeInjectionDiagnostics(const StartCells& start)
{
    return pImpl_->injDiag(start);
}

Opm::FlowDiagnosticsTool::Reverse
Opm::FlowDiagnosticsTool::
computeProductionDiagnostics(const StartCells& start)
{
    return pImpl_->prodDiag(start);
}
