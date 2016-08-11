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

#include <opm/flowdiagnostics/FlowDiagnosticsTool.hpp>
#include <opm/flowdiagnostics/reorder/tarjan.h>
#include <opm/flowdiagnostics/utility/AssembledConnections.hpp>
#include <cassert>
#include <iostream>

namespace Opm
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
            return std::vector<double>(pv_.size(), 0.0);
        }

        struct LocalSolution {
            CellSetValues tof;
            CellSetValues concentration;
        };

        LocalSolution solveLocal(const CellSet& startset)
        {
            static_cast<void>(startset);
            static_cast<void>(g_);
            static_cast<void>(flux_);
            return LocalSolution{ CellSetValues{}, CellSetValues{} };
        }

    private:
        void computeOrdering()
        {
            // Create the data structure needed for Tarjan's algorithm.
            const size_t num_cells = g_.numCells();
            const size_t num_connections = g_.numConnections();
            Utility::AssembledConnections ac;
            for (size_t conn_idx = 0; conn_idx < num_connections; ++conn_idx) {
                auto cells = g_.connection(conn_idx);
                const double connection_flux = flux_(ConnectionValues::ConnID{conn_idx},
                                                     ConnectionValues::PhaseID{0});
                if (connection_flux > 0.0) {
                    ac.addConnection(cells.first, cells.second, connection_flux);
                } else {
                    ac.addConnection(cells.second, cells.first, -connection_flux);
                }
            }
            ac.compress();

            // We might have to pad the start pointers if the last
            // cell(s) did not have outgoing fluxes, to get the
            // traditional format expected by tarjan().
            auto sp = ac.startPointers();
            if (sp.size() != num_cells + 1) {
                assert(sp.size() < num_cells + 1);
                sp.insert(sp.end(), num_cells + 1 - sp.size(), sp.back());
            }
            assert(sp.size() == num_cells + 1);

            // Compute topological ordering.
            struct Deleter { void operator()(TarjanSCCResult* x) { destroy_tarjan_sccresult(x); } };
            std::unique_ptr<TarjanSCCResult, Deleter> result(tarjan(num_cells, sp.data(), ac.neighbourhood().data()));

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
            static_cast<void>(cell);
#if 0
            // Compute flux terms.
            // Sources have zero tof, and therefore do not contribute
            // to upwind_term. Sinks on the other hand, must be added
            // to the downwind_flux (note sign change resulting from
            // different sign conventions: pos. source is injection,
            // pos. flux is outflow).
            double upwind_term = 0.0;
            double downwind_flux = std::max(-source_[cell], 0.0);
            for (int i = grid_.cell_facepos[cell]; i < grid_.cell_facepos[cell+1]; ++i) {
                int f = grid_.cell_faces[i];
                double flux;
                int other;
                // Compute cell flux
                if (cell == grid_.face_cells[2*f]) {
                    flux  = darcyflux_[f];
                    other = grid_.face_cells[2*f+1];
                } else {
                    flux  =-darcyflux_[f];
                    other = grid_.face_cells[2*f];
                }
                // Add flux to upwind_term or downwind_flux
                if (flux < 0.0) {
                    // Using tof == 0 on inflow, so we only add a
                    // nonzero contribution if we are on an internal
                    // face.
                    if (other != -1) {
                        upwind_term += flux*tof_[other];
                    }
                } else {
                    downwind_flux += flux;
                }
            }

            // Compute tof.
            tof_[cell] = (porevolume_[cell] - upwind_term)/downwind_flux;
#endif
        }


        void solveMultiCell(const int size, const int* cells)
        {
            static_cast<void>(size);
            static_cast<void>(cells);
        }

        const ConnectivityGraph& g_;
        const std::vector<double>& pv_;
        const ConnectionValues& flux_;
        std::vector<int> sequence_;
        std::vector<int> component_starts_;
    };

} // namespace Opm

#endif // OPM_TRACERTOFSOLVER_HEADER_INCLUDED
