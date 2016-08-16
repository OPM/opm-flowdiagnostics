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
#include <stdexcept>

namespace Opm
{
namespace FlowDiagnostics
{


    /// Class for solving the tracer and time-of-flight equations.
    ///
    /// Implements a first-order finite volume solver for
    /// (single-phase) time-of-flight using reordering.
    /// The equation solved is:
    ///     \f[v \cdot \nabla\tau = \phi\f]
    /// in which \f$ v \f$ is the fluid velocity, \f$ \tau \f$ is time-of-flight and
    /// \f$ \phi \f$ is the porosity. This is a boundary value problem, and
    /// \f$ \tau \f$ is considered zero on all inflow boundaries (or well inflows).
    ///
    /// The tracer equation is the same, except for the right hand side which is zero
    /// instead of \f[\phi\f].
    class TracerTofSolver
    {
    public:
        /// Initialize solver with a given flow graph (a weighted,
        /// directed asyclic graph) containing the out-fluxes from
        /// each cell, and pore volumes.
        TracerTofSolver(const AssembledConnections& graph,
                        const std::vector<double>& pore_volumes)
            : g_(graph),
              pv_(pore_volumes)
        {
        }



        /// Compute the global (combining all sources) time-of-flight of each cell.
        ///
        /// TODO: also compute tracer solution.
        std::vector<double> solveGlobal(const std::vector<CellSet>& all_startsets)
        {
            static_cast<void>(all_startsets); // TODO: use to compute tracer.

            // Reset instance variables.
            const int num_cells = pv_.size();
            upwind_influx_.clear();
            upwind_influx_.resize(num_cells, 0.0);
            upwind_contrib_.clear();
            upwind_contrib_.resize(num_cells, 0.0);
            tof_.clear();
            tof_.resize(num_cells, -1e100);
            num_multicell_ = 0;
            max_size_multicell_ = 0;
            max_iter_multicell_ = 0;

            // Compute topological ordering.
            computeOrdering();

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

            // Return computed time-of-flight.
            return tof_;
        }



        /// Output data struct for solveLocal().
        struct LocalSolution {
            CellSetValues tof;
            CellSetValues concentration;
        };



        /// Compute a local solution tracer and time-of-flight solution.
        ///
        /// Local means that only cells downwind from he startset are considered.
        /// The solution is therefore potentially sparse.
        /// TODO: not implemented!
        LocalSolution solveLocal(const CellSet& startset)
        {
            static_cast<void>(startset);
            return LocalSolution{ CellSetValues{}, CellSetValues{} };
        }



    private:

        // --------------  Private data members --------------

        const AssembledConnections& g_;
        const std::vector<double>& pv_;
        std::vector<int> sequence_;
        std::vector<int> component_starts_;
        std::vector<double> upwind_influx_;
        std::vector<double> upwind_contrib_;
        std::vector<double> tof_;
        int num_multicell_ = 0;
        int max_size_multicell_ = 0;
        int max_iter_multicell_ = 0;
        const double gauss_seidel_tol_ = 1e-3;

        // --------------  Private methods --------------

        void computeOrdering()
        {
            // We might have to pad the start pointers if the last
            // cell(s) did not have outgoing fluxes, to get the
            // traditional format expected by tarjan().
            const size_t num_cells = pv_.size();
            auto sp = g_.startPointers();
            if (sp.size() != num_cells + 1) {
                assert(sp.size() < num_cells + 1);
                sp.insert(sp.end(), num_cells + 1 - sp.size(), sp.back());
            }
            assert(sp.size() == num_cells + 1);

            // Compute reverse topological ordering.
            struct Deleter { void operator()(TarjanSCCResult* x) { destroy_tarjan_sccresult(x); } };
            std::unique_ptr<TarjanSCCResult, Deleter>
                result(tarjan(num_cells, sp.data(), g_.neighbourhood().data()));

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



        void solveSingleCell(const int cell)
        {
            // Compute downwind fluxes.
            double downwind_flux = 0.0;
            const auto& sp = g_.startPointers();
            if (cell < int(sp.size()) - 1) { // TODO: remove test when CRS from AC valid.
                for (const auto& conn : g_.cellNeighbourhood(cell)) {
                    downwind_flux += conn.weight;
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
                for (const auto& conn : g_.cellNeighbourhood(cell)) {
                    const int downwind_cell = conn.neighbour;
                    const double flux = conn.weight;
                    upwind_influx_[downwind_cell] += flux;
                    upwind_contrib_[downwind_cell] += tof_[cell] * flux;
                }
            }
        }



        void solveMultiCell(const int num_cells, const int* cells)
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
    };

} // namespace FlowDiagnostics
} // namespace Opm

#endif // OPM_TRACERTOFSOLVER_HEADER_INCLUDED
