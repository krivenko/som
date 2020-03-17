/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#include <triqs/utility/exceptions.hpp>

#include "observables.hpp"

namespace som {

using namespace triqs::arrays;
using namespace triqs::gfs;

bool is_stat_relevant(observable_kind kind) { return kind != ZeroTemp; }

statistic_enum observable_statistics(observable_kind kind) {
  switch(kind) {
    case FermionGf: return Fermion;
    case BosonCorr: return Boson;
    case BosonAutoCorr: return Boson;
    default:
      TRIQS_RUNTIME_ERROR << "observable_statistics: Statistics is undefined "
                             "for this observable kind";
  }
}

std::string observable_name(observable_kind kind) {
  switch(kind) {
    case FermionGf: return "Green's functions";
    case BosonCorr: return "bosonic correlator";
    case BosonAutoCorr: return "bosonic autocorrelator";
    case ZeroTemp: return "correlator at zero temperature";
    default: TRIQS_RUNTIME_ERROR << "observable_name: Invalid observable kind";
  }
}

std::pair<double, double> max_energy_window(observable_kind kind) {
  switch(kind) {
    case FermionGf: return std::make_pair(-HUGE_VAL, HUGE_VAL);
    case BosonCorr: return std::make_pair(-HUGE_VAL, HUGE_VAL);
    case BosonAutoCorr: return std::make_pair(0, HUGE_VAL);
    case ZeroTemp: return std::make_pair(0, HUGE_VAL);
    default:
      TRIQS_RUNTIME_ERROR << "max_energy_window: Invalid observable kind";
  }
}

// Construct a real-frequency GF from a configuration
void back_transform(observable_kind kind, configuration const& conf,
                    cache_index& ci, gf_view<refreq, scalar_valued> g_w) {
  bool bosoncorr = kind == BosonCorr || kind == BosonAutoCorr;

  g_w() = 0;

  for(auto const& rect : conf) {
    for(auto e : g_w.mesh())
      g_w.data()(e.index()) += rect.hilbert_transform(double(e), bosoncorr);

    if(kind == BosonAutoCorr) {
      // Add a reflected rectangle
      rectangle reflected_rect(-rect.center, rect.width, rect.height, ci);
      for(auto e : g_w.mesh())
        g_w.data()(e.index()) +=
            reflected_rect.hilbert_transform(double(e), true);
    }
  }

  if(bosoncorr) { g_w.data() *= -1.0 / M_PI; }
}

// Compute the GF tail from a configuration
triqs::arrays::array<std::complex<double>, 1>
compute_tail(observable_kind kind, configuration const& conf, cache_index& ci,
             int max_order) {
  assert(max_order >= 0);

  bool bosoncorr = kind == BosonCorr || kind == BosonAutoCorr;

  array<std::complex<double>, 1> tail(max_order + 1);
  tail() = .0;

  for(auto const& rect : conf) {
    tail += rect.tail_coefficients(0, max_order, bosoncorr);

    if(kind == BosonAutoCorr) {
      // Add a reflected rectangle
      rectangle reflected_rect(-rect.center, rect.width, rect.height, ci);
      tail += reflected_rect.tail_coefficients(0, max_order, true);
    }
  }

  if(bosoncorr) tail *= -1.0 / M_PI;

  return tail;
}

} // namespace som
