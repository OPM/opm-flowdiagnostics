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

#include <opm/flowdiagnostics/Toolbox.hpp>

#include <opm/flowdiagnostics/CellSet.hpp>
#include <opm/flowdiagnostics/ConnectionValues.hpp>
#include <opm/flowdiagnostics/ConnectivityGraph.hpp>
#include <opm/flowdiagnostics/TracerTofSolver.hpp>

#include <algorithm>
#include <exception>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <opm/utility/numeric/RandomVector.hpp>


namespace Opm
{
namespace FlowDiagnostics
{


// ---------------------------------------------------------------------
// Class Solution::Impl
// ---------------------------------------------------------------------

class Solution::Impl
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
Solution::Impl::assignToF(GlobalToF&& tof)
{
    tof_ = std::move(tof);
}

void
Solution::
Impl::assign(const CellSetID& i, TimeOfFlight&& tof)
{
    assign(i, std::move(tof.data), tracerToF_);
}

void
Solution::
Impl::assign(const CellSetID& i, Concentration&& conc)
{
    assign(i, std::move(conc.data), tracer_);
}

std::vector<CellSetID>
Solution::Impl::startPoints() const
{
    auto s = std::vector<CellSetID>{};
    s.reserve(tracer_.size());

    for (const auto& t : tracer_) {
        s.emplace_back(t.first);
    }

    return s;
}

const Solution::Impl::GlobalToF&
Solution::Impl::timeOfFlight() const
{
    return tof_;
}

CellSetValues
Solution::
Impl::timeOfFlight(const CellSetID& tracer) const
{
    return solutionValues(tracer, tracerToF_);
}

CellSetValues
Solution::
Impl::concentration(const CellSetID& tracer) const
{
    return solutionValues(tracer, tracer_);
}

void
Solution::
Impl::assign(const CellSetID& i,
             CellSetValues&&  x,
             SolutionMap&     soln)
{
    soln[i] = std::move(x);
}

CellSetValues
Solution::
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
// Class Toolbox::Impl
// ---------------------------------------------------------------------

class Toolbox::Impl
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

    AssembledConnections inj_conn_;
    AssembledConnections prod_conn_;
    bool conn_built_ = false;

    void buildAssembledConnections();
};

Toolbox::Impl::Impl(ConnectivityGraph g)
    : g_   (std::move(g))
    , pvol_(g_.numCells(), 0.0)
    , flux_(ConnectionValues::NumConnections{ 0 },
            ConnectionValues::NumPhases     { 0 })
{}

void
Toolbox::Impl::assign(const PoreVolume& pv)
{
    if (pv.data.size() != pvol_.size()) {
        throw std::logic_error("Inconsistently sized input "
                               "pore-volume field");
    }

    pvol_ = pv.data;
}

void
Toolbox::Impl::assign(const ConnectionFlux& flux)
{
    if (flux.data.numConnections() != g_.numConnections()) {
        throw std::logic_error("Inconsistently sized input "
                               "flux field");
    }

    flux_ = flux.data;
    conn_built_ = false;
}

Toolbox::Forward
Toolbox::Impl::injDiag(const StartCells& start)
{
    if (!conn_built_) {
        buildAssembledConnections();
    }

    using Soln       = Solution::Impl;
    using ToF        = Soln::TimeOfFlight;
    using Conc       = Soln::Concentration;

    using SolnPtr = std::unique_ptr<Soln>;

    SolnPtr x(new Soln());

    TracerTofSolver solver(inj_conn_, pvol_);
    x->assignToF(solver.solveGlobal(start.points));

    for (const auto& pt : start.points) {
        auto solution = solver.solveLocal(pt);
        x->assign(pt.id(), ToF{ solution.tof });
        x->assign(pt.id(), Conc{ solution.concentration });
    }

    return Forward{ Solution(std::move(x)) };
}

Toolbox::Reverse
Toolbox::Impl::prodDiag(const StartCells& start)
{
    if (!conn_built_) {
        buildAssembledConnections();
    }

    using Soln       = Solution::Impl;
    using ToF        = Soln::TimeOfFlight;
    using Conc       = Soln::Concentration;

    using SolnPtr = std::unique_ptr<Soln>;

    SolnPtr x(new Soln());

    TracerTofSolver solver(prod_conn_, pvol_);
    x->assignToF(solver.solveGlobal(start.points));

    for (const auto& pt : start.points) {
        auto solution = solver.solveLocal(pt);
        x->assign(pt.id(), ToF{ solution.tof });
        x->assign(pt.id(), Conc{ solution.concentration });
    }

    return Reverse{ Solution(std::move(x)) };
}

void
Toolbox::Impl::buildAssembledConnections()
{
    // Create the data structures needed by the tracer/tof solver.
    const size_t num_connections = g_.numConnections();
    const size_t num_phases = flux_.numPhases();
    inj_conn_ = AssembledConnections();
    prod_conn_ = AssembledConnections();
    for (size_t conn_idx = 0; conn_idx < num_connections; ++conn_idx) {
        auto cells = g_.connection(conn_idx);
        using ConnID = ConnectionValues::ConnID;
        using PhaseID = ConnectionValues::PhaseID;
        // Adding up all phase fluxes. TODO: ensure rigor, allow phase-based calculations.
        double connection_flux = 0.0;
        for (size_t phase = 0; phase < num_phases; ++phase) {
            connection_flux += flux_(ConnID{conn_idx}, PhaseID{phase});
        }
        if (connection_flux > 0.0) {
            inj_conn_.addConnection(cells.first, cells.second, connection_flux);
            prod_conn_.addConnection(cells.second, cells.first, connection_flux);
        } else {
            inj_conn_.addConnection(cells.second, cells.first, -connection_flux);
            prod_conn_.addConnection(cells.first, cells.second, -connection_flux);
        }
    }
    const int num_cells = g_.numCells();
    inj_conn_.compress(num_cells);
    prod_conn_.compress(num_cells);

    // Mark as built (until flux changed).
    conn_built_ = true;
}

// =====================================================================
// Implementation of public interface below separator
// =====================================================================

// ---------------------------------------------------------------------
// Class Solution
// ---------------------------------------------------------------------

Solution::~Solution()
{}

Solution::
Solution(const FlowDiagnostics::Solution& rhs)
    : pImpl_(new Impl(*rhs.pImpl_))
{}

Solution::
Solution(FlowDiagnostics::Solution&& rhs)
    : pImpl_(std::move(rhs.pImpl_))
{}

Solution::
Solution(std::unique_ptr<Impl> pImpl)
    : pImpl_(std::move(pImpl))
{}

std::vector<CellSetID>
Solution::startPoints() const
{
    return pImpl_->startPoints();
}

const std::vector<double>&
Solution::timeOfFlight() const
{
    return pImpl_->timeOfFlight();
}

CellSetValues
Solution::timeOfFlight(const CellSetID& tracer) const
{
    return pImpl_->timeOfFlight(tracer);
}

CellSetValues
Solution::concentration(const CellSetID& tracer) const
{
    return pImpl_->concentration(tracer);
}

// ---------------------------------------------------------------------
// Class Toolbox
// ---------------------------------------------------------------------

Toolbox::
Toolbox(const ConnectivityGraph& conn)
    : pImpl_(new Impl(conn))
{}

Toolbox::~Toolbox()
{}

Toolbox::Toolbox(Toolbox&& rhs)
    : pImpl_(std::move(rhs.pImpl_))
{
}

Toolbox&
Toolbox::operator=(Toolbox&& rhs)
{
    pImpl_ = std::move(rhs.pImpl_);

    return *this;
}

Toolbox&
Toolbox::assign(const PoreVolume& pv)
{
    pImpl_->assign(pv);

    return *this;
}

Toolbox&
Toolbox::assign(const ConnectionFlux& flux)
{
    pImpl_->assign(flux);

    return *this;
}

Toolbox::Forward
Toolbox::
computeInjectionDiagnostics(const StartCells& start)
{
    return pImpl_->injDiag(start);
}

Toolbox::Reverse
Toolbox::
computeProductionDiagnostics(const StartCells& start)
{
    return pImpl_->prodDiag(start);
}


} // namespace FlowDiagnostics
} // namespace Opm
