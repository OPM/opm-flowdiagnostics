/*
  Copyright 2015, 2016 SINTEF ICT, Applied Mathematics.

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

#include <opm/flowdiagnostics/DerivedQuantities.hpp>
#include <numeric>

namespace Opm
{
namespace FlowDiagnostics
{

    /// The F-Phi curve.
    ///
    /// The F-Phi curve is an analogue to the fractional flow
    /// curve in a 1D displacement. It can be used to compute
    /// other interesting diagnostic quantities such as the Lorenz
    /// coefficient. For a technical description see Shavali et
    /// al. (SPE 146446), Shook and Mitchell (SPE 124625).
    Graph flowCapacityStorageCapacityCurve(const Toolbox::Forward& injector_solution,
                                           const Toolbox::Reverse& producer_solution,
                                           const std::vector<double>& pv)
    {
        const auto& ftof = injector_solution.fd.timeOfFlight();
        const auto& rtof = producer_solution.fd.timeOfFlight();
        if (pv.size() != ftof.size() || pv.size() != rtof.size()) {
            throw std::runtime_error("computeFandPhi(): Input vectors must have same size.");
        }

        // Sort according to total travel time.
        const int n = pv.size();
        typedef std::pair<double, double> D2;
        std::vector<D2> time_and_pv(n);
        for (int ii = 0; ii < n; ++ii) {
            time_and_pv[ii].first = ftof[ii] + rtof[ii]; // Total travel time.
            time_and_pv[ii].second = pv[ii];
        }
        std::sort(time_and_pv.begin(), time_and_pv.end());

        // Compute Phi.
        std::vector<double> Phi(n + 1);
        Phi[0] = 0.0;
        for (int ii = 0; ii < n; ++ii) {
            Phi[ii+1] = time_and_pv[ii].second;
        }
        std::partial_sum(Phi.begin(), Phi.end(), Phi.begin());
        const double vt = Phi.back(); // Total pore volume.
        for (int ii = 1; ii < n+1; ++ii) { // Note limits of loop.
            Phi[ii] /= vt; // Normalize Phi.
        }

        // Compute F.
        std::vector<double> F(n + 1);
        F[0] = 0.0;
        for (int ii = 0; ii < n; ++ii) {
            F[ii+1] = time_and_pv[ii].second / time_and_pv[ii].first;
        }
        std::partial_sum(F.begin(), F.end(), F.begin());
        const double ft = F.back(); // Total flux.
        for (int ii = 1; ii < n+1; ++ii) { // Note limits of loop.
            F[ii] /= ft; // Normalize Phi.
        }

        return std::make_pair(F, Phi);
    }




    /// The Lorenz coefficient from the F-Phi curve.
    ///
    /// The Lorenz coefficient is a measure of heterogeneity. It
    /// is equal to twice the area between the F-Phi curve and the
    /// F = Phi line.  The coefficient can vary from zero to
    /// one. If the coefficient is zero (so the F-Phi curve is a
    /// straight line) we have perfect piston-like displacement
    /// while a coefficient of one indicates infinitely
    /// heterogenous displacement (essentially no sweep).
    ///
    /// Note: The coefficient is analogous to the Gini coefficient
    /// of economic theory, where the name Lorenz curve is applied
    /// to what we call the F-Phi curve.
    double lorenzCoefficient(const Graph& )//flowcap_storagecap_curve)
    {
        return double();
    }




    /// Compute sweep efficiency versus dimensionless time (pore
    /// volumes injected).
    ///
    /// The sweep efficiency is analogue to 1D displacement using
    /// the F-Phi curve as flux function.
    Graph sweepEfficiency(const Toolbox::Forward& ,//injector_solution,
                          const Toolbox::Reverse& ,//producer_solution,
                          const std::vector<double>& )//pore_volume)
    {
        return Graph();
    }



} // namespace FlowDiagnostics
} // namespace Opm
