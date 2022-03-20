/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * SOM is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * SOM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * SOM. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <utility>
#include <vector>

#include <nda/nda.hpp>
#include <triqs/mesh/refreq.hpp>

#include "som_core/som_core.hpp"

// Functions in this header implement a statistical analysis technique for
// spectral functions described in Sections I-II of
//
// [1] O. Goulko et al. Phys. Rev. B 95, 014102 (2017).

namespace som {

// Type of resolution function \bar K(m, z)
enum resolution_function : unsigned int {
  rectangle,  // Rectangle R(z_m, \Delta_m, 1/\Delta_m; z)
  lorentzian, // \frac{1}{\pi} \frac{\Delta_m/2}{(z - z_m)^2 + (\Delta_m/2)^2}
  gaussian    // \frac{1}{\sqrt{2\pi(\Delta_m/2)^2}}
              //    \exp(-\frac{(z - z_m)^2}{2(\Delta_m/2)^2})
};

/// Spectral integral i_m^{(j)} (Eq. (4) of [1])
///
/// i_m^{(j)} = \int_{-\infty}^\infty dz \bar K(m, z) A^{(j)}(z)
double spectral_integral(double z_m,
                         double delta_m,
                         configuration const& c,
                         resolution_function r_func);

/// Spectral integrals i_m^{(j)} (Eq. (4) of [1]) for all points z_m
/// of a regular real frequency mesh.
nda::vector<double> spectral_integral(triqs::mesh::refreq const& mesh,
                                      configuration const& c,
                                      resolution_function r_func);

/// Spectral integrals i_m^{(j)} (Eq. (4) of [1]) for an arbitrary set
/// of real frequency intervals.
nda::vector<double>
spectral_integral(std::vector<std::pair<double, double>> const& intervals,
                  configuration const& c,
                  resolution_function r_func);

/// Spectral averages i_m (Eq. (6) of [1])
///
/// i_m = J^{-1} \sum_{j=1}^J i_m^{(j)}

// Regular real frequency mesh
nda::vector<double> spectral_avg(som_core const& cont,
                                 long i,
                                 triqs::mesh::refreq const& mesh,
                                 resolution_function r_func);

// Arbitrary set of real frequency intervals
nda::vector<double>
spectral_avg(som_core const& cont,
             long i,
             std::vector<std::pair<double, double>> const& intervals,
             resolution_function r_func);

/// Spectral dispersions \sigma_m (Eq. (7) of [1])
///
/// \sigma_m^2 = J^{-1} \usm_{j=1}^J (i_m^{(j)} - i_m)^2

// Regular real frequency mesh
nda::vector<double> spectral_disp(som_core const& cont,
                                  long i,
                                  triqs::gfs::refreq const& mesh,
                                  nda::vector<double> const& avg,
                                  resolution_function r_func);

// Arbitrary set of real frequency intervals
nda::vector<double>
spectral_disp(som_core const& cont,
              long i,
              std::vector<std::pair<double, double>> const& intervals,
              nda::vector<double> const& avg,
              resolution_function r_func);

/// Spectral correlation matrix \sigma_{mm'} (Eq. (8) of [1])
///
/// \sigma_{mm'} = J^{-1} \usm_{j=1}^J (i_m^{(j)} - i_m)(i_{m'}^{(j)} - i_{m'})

// Regular real frequency mesh
nda::matrix<double> spectral_corr(som_core const& cont,
                                  long i,
                                  triqs::mesh::refreq const& mesh,
                                  nda::vector<double> const& avg,
                                  resolution_function r_func);

// Arbitrary set of real frequency intervals
nda::matrix<double>
spectral_corr(som_core const& cont,
              long i,
              std::vector<std::pair<double, double>> const& intervals,
              nda::vector<double> const& avg,
              resolution_function r_func);

} // namespace som
