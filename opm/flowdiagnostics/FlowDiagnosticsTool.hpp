/*
  Copyright 2016 Statoil ASA.
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

#ifndef OPM_FLOWDIAGNOSTICSTOOL_HEADER_INCLUDED
#define OPM_FLOWDIAGNOSTICSTOOL_HEADER_INCLUDED

#include <opm/flowdiagnostics/CellSet.hpp>
#include <opm/flowdiagnostics/CellSetValues.hpp>

#include <memory>
#include <vector>

namespace Opm
{
    // These classes are not defined for the moment, definitions will
    // depend on needs and be refined as we progress.

    // Nodes are cells, edges are connections.  Wells not included.
    class ConnectivityGraph;

    // One element per phase for each connection.
    class ConnectionValues;

    /// Toolbox for running flow diagnostics.
    class FlowDiagnosticsTool;

    /// Results from diagnostics computations.
    class FlowDiagnosticsSolution
    {
    public:
        ~FlowDiagnosticsSolution();

        FlowDiagnosticsSolution(const FlowDiagnosticsSolution& rhs);
        FlowDiagnosticsSolution(FlowDiagnosticsSolution&& rhs);

        /// Ids of stored tracer solutions.
        std::vector<CellSetID> startPoints() const;

        /// Time-of-flight field from all start points.
        const std::vector<double>& timeOfFlight() const;

        /// Time-of-flight field restricted to single tracer region.
        CellSetValues timeOfFlight(const CellSetID& tracer) const;

        /// The computed tracer field corresponding to a single tracer.
        ///
        /// The \c tracer must correspond to an id passed in
        /// computeX...Diagnostics().
        CellSetValues concentration(const CellSetID& tracer) const;

        friend class FlowDiagnosticsTool;

    private:
        class Impl;

        explicit FlowDiagnosticsSolution(std::unique_ptr<Impl> pImpl);

        std::unique_ptr<Impl> pImpl_;
    };

    /// Toolbox for running flow diagnostics.
    class FlowDiagnosticsTool
    {
    public:
        /// Construct from known neighbourship relation.
        explicit FlowDiagnosticsTool(const ConnectivityGraph& connectivity);

        ~FlowDiagnosticsTool();

        struct PoreVolume
        {
            const std::vector<double>& data;
        };

        struct ConnectionFlux
        {
            const ConnectionValues& data;
        };

        struct StartCells
        {
            const std::vector<CellSet>& points;
        };

        struct Forward
        {
            const FlowDiagnosticsSolution fd;
        };

        struct Reverse
        {
            const FlowDiagnosticsSolution fd;
        };

        FlowDiagnosticsTool& assign(const PoreVolume&     pv);
        FlowDiagnosticsTool& assign(const ConnectionFlux& flux);

        /// Compute forward time-of-flight and tracer solutions.
        ///
        /// An element of \code start.points \endcode provides a set of
        /// starting locations for a single tracer.
        ///
        /// Forward time-of-flight is the time needed for a neutral fluid
        /// particle to flow from the nearest fluid source to an arbitrary
        /// point in the model.  The tracer solutions identify cells that
        /// are flooded by injectors.
        ///
        /// The IDs of the \code start.points \endcode must be unique.
        Forward computeInjectionDiagnostics(const StartCells& start);

        /// Compute reverse time-of-flight and tracer solutions.
        ///
        /// An element of \code start.points \endcode provides a set of
        /// starting locations for a single tracer.
        ///
        /// Reverse time-of-flight is the time needed for a neutral fluid
        /// particle to flow from an arbitrary point to the nearest fluid
        /// sink in the model.  The tracer solutions identify cells that are
        /// drained by producers.
        ///
        /// The IDs of the \code start.points \endcode must be unique.
        Reverse computeProductionDiagnostics(const StartCells& start);

    private:
        class Impl;

        std::unique_ptr<Impl> pImpl_;
    };
} // namespace Opm

#endif // OPM_FLOWDIAGNOSTICSTOOL_HEADER_INCLUDED
