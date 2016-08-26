/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.

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

#if HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <opm/flowdiagnostics/TracerTofSolver.hpp>
#include <opm/utility/graph/tarjan.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

namespace Opm
{
namespace FlowDiagnostics
{




    // This computes both in and outflux with a single traversal of the graph.
    struct TracerTofSolver::InOutFluxComputer
    {
        InOutFluxComputer(const AssembledConnections& graph)
        {
            const int num_cells = graph.numRows();
            influx.resize(num_cells, 0.0);
            outflux.resize(num_cells, 0.0);
            for (int cell = 0; cell < num_cells; ++cell) {
                const auto nb = graph.cellNeighbourhood(cell);
                for (const auto& conn : nb) {
                    influx[conn.neighbour] += conn.weight;
                    outflux[cell] += conn.weight;
                }
            }
        }

        std::vector<double> influx;
        std::vector<double> outflux;
    };





    TracerTofSolver::TracerTofSolver(const AssembledConnections& graph,
                                     const std::vector<double>& pore_volumes)
        : TracerTofSolver(graph, pore_volumes, InOutFluxComputer(graph))
    {
    }





    // The InOutFluxComputer is used so that influx_ and outflux_ can be
    // const members of the class.
    TracerTofSolver::TracerTofSolver(const AssembledConnections& graph,
                                     const std::vector<double>& pore_volumes,
                                     InOutFluxComputer&& inout)
        : g_(graph)
        , pv_(pore_volumes)
        , influx_(std::move(inout.influx))
        , outflux_(std::move(inout.outflux))
    {
    }





    std::vector<double> TracerTofSolver::solveGlobal(const std::vector<CellSet>& all_startsets)
    {
        // Reset solver variables and set source terms.
        prepareForSolve();
        for (const CellSet& startset : all_startsets) {
            setupSourceTerms(startset);
        }

        // Compute topological ordering and solve.
        computeOrdering();
        solve();

        // Return computed time-of-flight.
        return tof_;
    }





    TracerTofSolver::LocalSolution TracerTofSolver::solveLocal(const CellSet& startset)
    {
        // Reset solver variables and set source terms.
        prepareForSolve();
        setupSourceTerms(startset);

        // Compute local topological ordering and solve.
        computeLocalOrdering(startset);
        solve();

        // Return computed time-of-flight.
        CellSetValues local_tof;
        const int num_elements = component_starts_.back();
        for (int element = 0; element < num_elements; ++element) {
            const int cell = sequence_[element];
            local_tof.addCellValue(cell, tof_[cell]);
        }
        return LocalSolution{ local_tof, CellSetValues{} }; // TODO also return tracer
    }





    void TracerTofSolver::prepareForSolve()
    {
        // Reset instance variables.
        const int num_cells = pv_.size();
        source_term_.clear();
        source_term_.resize(num_cells, 0.0);
        upwind_contrib_.clear();
        upwind_contrib_.resize(num_cells, 0.0);
        tof_.clear();
        tof_.resize(num_cells, -1e100);
        num_multicell_ = 0;
        max_size_multicell_ = 0;
        max_iter_multicell_ = 0;
    }





    void TracerTofSolver::setupSourceTerms(const CellSet& startset)
    {
        for (const int cell : startset) {
            source_term_[cell] = std::numeric_limits<double>::infinity(); // std::max(outflux_[cell] - influx_[cell], 0.0); // Ignore negative sources.
        }
    }





    void TracerTofSolver::computeOrdering()
    {
        // Compute reverse topological ordering.
        const size_t num_cells = pv_.size();
        assert(g_.startPointers().size() == num_cells + 1);
        struct Deleter { void operator()(TarjanSCCResult* x) { destroy_tarjan_sccresult(x); } };
        std::unique_ptr<TarjanSCCResult, Deleter>
            result(tarjan(num_cells, g_.startPointers().data(), g_.neighbourhood().data()));

        // Must reverse ordering, since Tarjan computes reverse ordering.
        const int ok = tarjan_reverse_sccresult(result.get());
        if (!ok) {
            throw std::runtime_error("Failed to reverse topological ordering in TracerTofSolver::computeOrdering()");
        }

        // Extract data from solution.
        sequence_.resize(num_cells);
        const int num_comp = tarjan_get_numcomponents(result.get());
        component_starts_.resize(num_comp + 1);
        component_starts_[0] = 0;
        for (int comp = 0; comp < num_comp; ++comp) {
            const TarjanComponent tc = tarjan_get_strongcomponent(result.get(), comp);
            std::copy(tc.vertex, tc.vertex + tc.size, sequence_.begin() + component_starts_[comp]);
            component_starts_[comp + 1] = component_starts_[comp] + tc.size;
        }
        assert(component_starts_.back() == int(num_cells));
    }





    void TracerTofSolver::computeLocalOrdering(const CellSet& startset)
    {
        // Extract start cells.
        std::vector<int> startcells(startset.begin(), startset.end());

        // Compute reverse topological ordering.
        const size_t num_cells = pv_.size();
        assert(g_.startPointers().size() == num_cells + 1);
        struct ResultDeleter { void operator()(TarjanSCCResult* x) { destroy_tarjan_sccresult(x); } };
        std::unique_ptr<TarjanSCCResult, ResultDeleter> result;
        {
            struct WorkspaceDeleter { void operator()(TarjanWorkSpace* x) { destroy_tarjan_workspace(x); } };
            std::unique_ptr<TarjanWorkSpace, WorkspaceDeleter> ws(create_tarjan_workspace(num_cells));
            result.reset(tarjan_reachable_sccs(num_cells, g_.startPointers().data(), g_.neighbourhood().data(),
                                               startcells.size(), startcells.data(), ws.get()));
        }

        // Must reverse ordering, since Tarjan computes reverse ordering.
        const int ok = tarjan_reverse_sccresult(result.get());
        if (!ok) {
            throw std::runtime_error("Failed to reverse topological ordering in TracerTofSolver::computeOrdering()");
        }

        // Extract data from solution.
        sequence_.resize(num_cells);
        const int num_comp = tarjan_get_numcomponents(result.get());
        component_starts_.resize(num_comp + 1);
        component_starts_[0] = 0;
        for (int comp = 0; comp < num_comp; ++comp) {
            const TarjanComponent tc = tarjan_get_strongcomponent(result.get(), comp);
            std::copy(tc.vertex, tc.vertex + tc.size, sequence_.begin() + component_starts_[comp]);
            component_starts_[comp + 1] = component_starts_[comp] + tc.size;
        }
    }





    void TracerTofSolver::solve()
    {
        // Solve each component.
        const int num_components = component_starts_.size() - 1;
        for (int comp = 0; comp < num_components; ++comp) {
            const int comp_size = component_starts_[comp + 1] - component_starts_[comp];
            if (comp_size == 1) {
                solveSingleCell(sequence_[component_starts_[comp]]);
            } else {
                solveMultiCell(comp_size, &sequence_[component_starts_[comp]]);
            }
        }
    }





    void TracerTofSolver::solveSingleCell(const int cell)
    {
        // Note: comment below is not matching the current implementation,
        // but kept as a reference since we have not settled the formulation
        // 100% yet.
        //
        // If source cell, we have influx not accounted for, and
        // downwind_flux > upwind_flux_[cell]. However the tof at
        // the influx (wells typically) is zero, so there is no
        // contribution. However, we will customary halve the tof
        // for that cell.
        // If sink cell, we have outflux not accounted for, and
        // downwind_flux < upwind_flux_[cell]. In that case we
        // account for it by assuming incompressibility and making
        // them equal.
        //const double downwind_flux = std::max(outflux_[cell], influx_[cell]);

        const double total_influx_ = influx_[cell] + source_term_[cell];

        // Compute tof.
        tof_[cell] = (pv_[cell] + upwind_contrib_[cell])/total_influx_;

        // Halve if source (well inflow) cell.
        // if (upwind_contrib_[cell] == 0.0) {
        //     tof_[cell] = 0.5 * tof_[cell];
        // }

        // Set tof to zero if cell has no reservoir influx (assume it is a source).
        if (influx_[cell] == 0.0) {
            tof_[cell] = 0.0;
        }

        // Set contribution for my downwind cells (if any).
        for (const auto& conn : g_.cellNeighbourhood(cell)) {
            const int downwind_cell = conn.neighbour;
            const double flux = conn.weight;
            upwind_contrib_[downwind_cell] += tof_[cell] * flux;
        }
    }





    void TracerTofSolver::solveMultiCell(const int num_cells, const int* cells)
    {
        ++num_multicell_;
        max_size_multicell_ = std::max(max_size_multicell_, num_cells);
        // std::cout << "Multiblock solve with " << num_cells << " cells." << std::endl;

        // Using a Gauss-Seidel approach.
        double max_delta = 1e100;
        int num_iter = 0;
        while (max_delta > gauss_seidel_tol_) {
            max_delta = 0.0;
            ++num_iter;
            for (int ci = 0; ci < num_cells; ++ci) {
                const int cell = cells[ci];
                const double tof_before = tof_[cell];
                solveSingleCell(cell);
                max_delta = std::max(max_delta, std::fabs(tof_[cell] - tof_before));
            }
            // std::cout << "Max delta = " << max_delta << std::endl;
        }
        max_iter_multicell_ = std::max(max_iter_multicell_, num_iter);
    }




} // namespace FlowDiagnostics
} // namespace Opm
