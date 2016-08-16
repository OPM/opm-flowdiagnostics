/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_TRACERTOFSOLVER_HEADER_INCLUDED
#define OPM_TRACERTOFSOLVER_HEADER_INCLUDED

#include <opm/flowdiagnostics/Toolbox.hpp>
#include <opm/utility/graph/tarjan.h>
#include <opm/utility/graph/AssembledConnections.hpp>
#include <cassert>
#include <iostream>

namespace Opm
{
namespace FlowDiagnostics
{

    class TracerTofSolver
    {
    public:
        TracerTofSolver(const ConnectivityGraph& graph,
                        const std::vector<double>& pore_volumes,
                        const ConnectionValues& flux)
            : g_(graph),
              pv_(pore_volumes),
              flux_(flux)
        {
        }

        std::vector<double> solveGlobal(const std::vector<CellSet>& all_startsets)
        {
            static_cast<void>(all_startsets);
            upwind_influx_.clear();
            upwind_influx_.resize(g_.numCells(), 0.0);
            upwind_contrib_.clear();
            upwind_contrib_.resize(g_.numCells(), 0.0);
            tof_.clear();
            tof_.resize(g_.numCells(), -1e100);
            computeOrdering();
            const int num_components = component_starts_.size() - 1;
            for (int comp = 0; comp < num_components; ++comp) {
                const int comp_size = component_starts_[comp + 1] - component_starts_[comp];
                if (comp_size == 1) {
                    solveSingleCell(sequence_[component_starts_[comp]]);
                } else {
                    solveMultiCell(comp_size, &sequence_[component_starts_[comp]]);
                }
            }
            return tof_;
        }

        struct LocalSolution {
            CellSetValues tof;
            CellSetValues concentration;
        };

        LocalSolution solveLocal(const CellSet& startset)
        {
            static_cast<void>(startset);
            return LocalSolution{ CellSetValues{}, CellSetValues{} };
        }

    private:
        void computeOrdering()
        {
            // Create the data structure needed for Tarjan's algorithm.
            const size_t num_cells = g_.numCells();
            const size_t num_connections = g_.numConnections();
            assembled_conn_ = AssembledConnections();
            for (size_t conn_idx = 0; conn_idx < num_connections; ++conn_idx) {
                auto cells = g_.connection(conn_idx);
                const double connection_flux = flux_(ConnectionValues::ConnID{conn_idx},
                                                     ConnectionValues::PhaseID{0});
                if (connection_flux > 0.0) {
                    assembled_conn_.addConnection(cells.first, cells.second, connection_flux);
                } else {
                    assembled_conn_.addConnection(cells.second, cells.first, -connection_flux);
                }
            }
            assembled_conn_.compress();

            // We might have to pad the start pointers if the last
            // cell(s) did not have outgoing fluxes, to get the
            // traditional format expected by tarjan().
            auto sp = assembled_conn_.startPointers();
            if (sp.size() != num_cells + 1) {
                assert(sp.size() < num_cells + 1);
                sp.insert(sp.end(), num_cells + 1 - sp.size(), sp.back());
            }
            assert(sp.size() == num_cells + 1);

            // Compute reverse topological ordering.
            struct Deleter { void operator()(TarjanSCCResult* x) { destroy_tarjan_sccresult(x); } };
            std::unique_ptr<TarjanSCCResult, Deleter> result(tarjan(num_cells,
                                                                    sp.data(),
                                                                    assembled_conn_.neighbourhood().data()));

            // Must reverse ordering, since Tarjan computes reverse ordering.
            const int ok = tarjan_reverse_sccresult(result.get());
            if (!ok) {
                throw std::runtime_error("Failed to reverse ordering.");
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

        void solveSingleCell(const int cell)
        {
            // Compute downwind fluxes.
            double downwind_flux = 0.0;
            const auto& sp = assembled_conn_.startPointers();
            if (cell < int(sp.size()) - 1) { // TODO: remove test when CRS from AC valid.
                for (size_t connection_index = sp[cell]; connection_index < sp[cell + 1]; ++connection_index) {
                    const double flux = assembled_conn_.connectionWeight()[connection_index];
                    downwind_flux += flux;
                }
            }

            // If source cell, we have influx not accounted for, and
            // downwind_flux > upwind_flux_[cell]. However the tof at
            // the influx (wells typically) is zero, so there is no
            // contribution. However, we will customary halve the tof
            // for that cell.
            // If sink cell, we have outflux not accounted for, and
            // downwind_flux < upwind_flux_[cell]. In that case we
            // account for it by assuming incompressibility and making
            // them equal.
            downwind_flux = std::max(downwind_flux, upwind_influx_[cell]);

            // Compute tof.
            tof_[cell] = (pv_[cell] + upwind_contrib_[cell])/downwind_flux;

            // Halve if source (well inflow) cell.
            if (upwind_influx_[cell] == 0.0) {
                tof_[cell] = 0.5 * tof_[cell];
            }

            // Set contribution for my downwind cells (if any).
            if (cell < int(sp.size()) - 1) { // TODO: remove test when CRS from AC valid.
                for (size_t connection_index = sp[cell]; connection_index < sp[cell + 1]; ++connection_index) {
                    const int downwind_cell = assembled_conn_.neighbourhood()[connection_index];
                    const double flux = assembled_conn_.connectionWeight()[connection_index];
                    upwind_influx_[downwind_cell] += flux;
                    upwind_contrib_[downwind_cell] += tof_[cell] * flux;
                }
            }
        }


        void solveMultiCell(const int size, const int* cells)
        {
            static_cast<void>(size);
            static_cast<void>(cells);
        }

        const ConnectivityGraph& g_;
        const std::vector<double>& pv_;
        const ConnectionValues& flux_;
        AssembledConnections assembled_conn_;
        std::vector<int> sequence_;
        std::vector<int> component_starts_;
        std::vector<double> upwind_influx_;
        std::vector<double> upwind_contrib_;
        std::vector<double> tof_;
    };

} // namespace FlowDiagnostics
} // namespace Opm

#endif // OPM_TRACERTOFSOLVER_HEADER_INCLUDED
