/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2023 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#include <nda/nda.hpp>

#include "fit_quality.hpp"

namespace som {

template <typename KernelType>
fit_quality<KernelType>::fit_quality(KernelType const& kern,
                                     rhs_type const& rhs,
                                     rhs_type const& error_bars)
   : kern(kern), rhs(rhs), error_bars(error_bars) {}

template <typename KernelType>
double fit_quality<KernelType>::corr(double x1, double x2) {
  return (x1 * x2 < 0) ? 1 : 0;
}

template <typename KernelType>
double fit_quality<KernelType>::corr(std::complex<double> z1,
                                     std::complex<double> z2) {
  std::complex<double> w = z1 * std::conj(z2);
  return 0.5 * (1.0 - (w == .0 ? 1.0 : (w.real() / std::abs(w))));
}

template <typename KernelType>
double fit_quality<KernelType>::operator()(configuration const& c) const {
  rhs_type delta = (rhs - kern(c)) / error_bars;
  int M = delta.shape()[0];
  double kappa = 0;
  for(int i = 1; i < M; ++i) kappa += corr(delta(i), delta(i - 1));
  kappa /= M - 1;
  return kappa;
}

INSTANTIATE_CLASS_FOR_EACH_KERNEL(fit_quality)

} // namespace som
