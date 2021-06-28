/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2021 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#include <complex>
#include <string>
#include <utility>

#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/size.hpp>

#include <mpi/mpi.hpp>
#include <triqs/arrays.hpp>
#include <triqs/gfs.hpp>

#include "../config_update.hpp"
#include "../configuration.hpp"

namespace som {

// All supported kinds of observables
#define ALL_OBSERVABLES (FermionGf)(BosonCorr)(BosonAutoCorr)(ZeroTemp)

enum observable_kind : unsigned int { BOOST_PP_SEQ_ENUM(ALL_OBSERVABLES) };
constexpr const unsigned int n_observable_kinds =
    BOOST_PP_SEQ_SIZE(ALL_OBSERVABLES);

// Is statistics defined for this observable?
bool is_stat_relevant(observable_kind kind);

// Statistics of observables
triqs::gfs::statistic_enum observable_statistics(observable_kind kind);

// Names of observables
std::string observable_name(observable_kind kind);

// Widest energy windows for observables
std::pair<double, double> max_energy_window(observable_kind kind);

// Construct a real-frequency GF from a configuration
void back_transform(
    observable_kind kind, configuration const& conf,
    triqs::gfs::gf_view<triqs::gfs::refreq, triqs::gfs::scalar_valued> g_w,
    mpi::communicator const& comm = {});

// Compute the GF tail from a configuration
triqs::arrays::array<std::complex<double>, 1>
compute_tail(observable_kind kind, configuration const& conf, cache_index& ci,
             mpi::communicator& comm, int max_order);

} // namespace som
