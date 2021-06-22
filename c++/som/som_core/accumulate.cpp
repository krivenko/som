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

#include <cmath>
#include <limits>

#include <boost/preprocessor/seq/for_each_product.hpp>

#include <triqs/utility/signal_handler.hpp>

#include <som/kernels/all.hpp>
#include <som/objective_function.hpp>
#include <som/solution_worker.hpp>

#include "som_core.hpp"
#include "common.hxx"

namespace som {

using triqs::statistics::histogram;

/////////////////////
// som_core::run() //
/////////////////////

void validate_params(observable_kind kind, run_parameters_t & params) {
  double e_min, e_max;
  std::tie(e_min, e_max) = max_energy_window(kind);
  if(params.energy_window.first < e_min) {
    params.energy_window.first = e_min;
    if(params.verbosity > 0)
      warning("left boundary of the energy window is reset to " +
              std::to_string(e_min));
  }
  if(params.energy_window.second > e_max) {
    params.energy_window.second = e_max;
    if(params.verbosity > 0)
      warning("right boundary of the energy window is reset to " +
              std::to_string(e_max));
  }

  if(params.energy_window.first >= params.energy_window.second)
    fatal_error("wrong energy window [" +
                to_string(params.energy_window.first) + ";" +
                to_string(params.energy_window.second) + "]");

  if(params.min_rect_width <= 0 || params.min_rect_width > 1)
    fatal_error("min_rect_width must be in (0;1]");

  if(params.min_rect_weight <= 0 || params.min_rect_weight > 1)
    fatal_error("min_rect_weight must be in (0;1]");

  if(params.adjust_l_range.first > params.adjust_l_range.second)
    fatal_error("Wrong adjust_l_range");
}

void som_core::run(run_parameters_t const& p) {

  params = p;
  validate_params(kind, params);

  triqs::signal_handler::start();
  run_status = 0;
  try {
#define RUN_IMPL_CASE(r, okmk)                                                 \
  case(int(BOOST_PP_SEQ_ELEM(0, okmk)) +                                       \
       n_observable_kinds * mesh_traits<BOOST_PP_SEQ_ELEM(1, okmk)>::index):   \
    run_impl<kernel<BOOST_PP_SEQ_ENUM(okmk)>>();                               \
    break;
    switch(int(kind) + n_observable_kinds * mesh.index()) {
      BOOST_PP_SEQ_FOR_EACH_PRODUCT(RUN_IMPL_CASE,
                                    (ALL_OBSERVABLES)(ALL_INPUT_MESHES))
    }
#undef RUN_IMPL_CASE
  } catch(stopped& e) {
    run_status = e.code;
    triqs::signal_handler::received(true);
  }
  triqs::signal_handler::stop();
}

template <typename KernelType> void som_core::run_impl() {

  auto stop_callback = triqs::utility::clock_callback(params.max_time);

  using mesh_t = typename KernelType::mesh_type;
  mesh_t const& m = std::get<mesh_t>(mesh);

  if(params.verbosity > 0) {
    std::cout << "Constructing integral kernel... " << std::flush;
  }
  KernelType kernel(m);
  if(params.verbosity > 0) {
    std::cout << "done" << std::endl;
    std::cout << "Kernel: " << kernel << std::endl;
  }

  // Find solution for each component of GF
  for(int i = 0; i < data.size(); ++i) {
    auto & d = data[i];

    // Prepare histogram
    if(params.make_histograms)
      d.histogram = histogram();

    if(params.verbosity > 0)
      std::cout << "Running algorithm for observable component [" << i << ","
                << i << "]" << std::endl;

    d.final_solution = accumulate(kernel, d, stop_callback, params.f);
  }

  ci.invalidate_all();
}

template <typename KernelType>
configuration som_core::accumulate(KernelType const& kern,
                                   data_t & d,
                                   std::function<bool()> const& stop_callback,
                                   int F) {

  if(params.verbosity >= 1)
    std::cout << "Accumulating particular solutions ..." << std::endl;

  using mesh_t = typename KernelType::mesh_type;
  auto const& rhs = d.get_rhs<mesh_t>();
  auto const& error_bars = d.get_error_bars<mesh_t>();

  objective_function<KernelType> of(kern, rhs, error_bars);
  solution_worker<KernelType> worker(of,
                                     d.norm,
                                     ci,
                                     params,
                                     stop_callback,
                                     F);
  auto& rng = worker.get_rng();

  int n_sol_max = 0; // Number of solutions to be accumulated
  int n_sol, i = 0;  // Global and rank-local indices of solution
  int n_good_solutions,
      n_verygood_solutions;   // Number of good and very good solutions
  double objf_min = HUGE_VAL; // Minimal value of D
  do {
    if(params.adjust_l) {
      n_sol_max += params.adjust_l_range.first;
      if(n_sol_max > params.adjust_l_range.second) {
        if(params.verbosity >= 1)
          warning("Upper bound of adjust_l_range has been reached");
        break;
      }
    } else
      n_sol_max = params.l;
    d.basis_solutions.clear();
    d.basis_solutions.reserve(n_sol_max);

    if(params.verbosity >= 1)
      std::cout
          << "Increasing the total number of solutions to be accumulated to "
          << n_sol_max << std::endl;

    for(; (n_sol = comm.rank() + i * comm.size()) < n_sol_max; ++i) {
      if(params.verbosity >= 2) {
        std::cout << "[Rank " << comm.rank()
                  << "] Accumulation of particular solution " << n_sol
                  << std::endl;
      }

      d.basis_solutions.emplace_back(worker(1 + rng(params.max_rects)), 0);

      double D = worker.get_objf_value();
      d.basis_solutions.back().second = D;
      objf_min = std::min(objf_min, D);

      if(params.verbosity >= 2) {
        std::cout << "[Rank " << comm.rank() << "] Solution " << n_sol
                  << ": D = " << D << std::endl;
      }
    }
    comm.barrier();

    // Global minimum of D_min
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    objf_min = mpi::all_reduce(objf_min, comm, 0, MPI_MIN);

    // Recalculate numbers of good and very good solutions
    n_good_solutions = n_verygood_solutions = 0;
    for(auto const& s : d.basis_solutions) {
      if(s.second / objf_min <= params.adjust_l_good_d) ++n_good_solutions;
      if(s.second / objf_min <= params.adjust_l_verygood_d)
        ++n_verygood_solutions;
    }
    n_good_solutions = mpi::all_reduce(n_good_solutions);
    n_verygood_solutions = mpi::all_reduce(n_verygood_solutions);

    if(params.verbosity >= 1) {
      std::cout << "D_min = " << objf_min << std::endl;
      std::cout << "Number of good solutions (D/D_min <= "
                << params.adjust_l_good_d << ") = " << n_good_solutions
                << std::endl;
      std::cout << "Number of very good solutions (D/D_min <= "
                << params.adjust_l_verygood_d << ") = " << n_verygood_solutions
                << std::endl;
    }

  } while(params.adjust_l &&
          double(n_verygood_solutions) / double(n_good_solutions) <
              params.adjust_l_ratio);

  comm.barrier();

  if(params.verbosity >= 1) {
    std::cout << "Accumulation complete." << std::endl;
    std::cout << "Summing up good solutions ..." << std::endl;
  }

  if(d.histogram)
    *d.histogram = histogram(objf_min,
                             objf_min * params.hist_max,
                             params.hist_n_bins);

  configuration sol_sum(ci);

  // Rank-local stage of summation
  for(auto const& s : d.basis_solutions) {
    if(d.histogram) *d.histogram << s.second;
    // Pick only good solutions
    if(s.second / objf_min <= params.adjust_l_good_d) sol_sum += s.first;
  }
  sol_sum *= 1.0 / double(n_good_solutions);

  // Sum over all processes
  sol_sum = mpi_reduce(sol_sum, comm, 0, true);

  if(d.histogram)
    *d.histogram = mpi_reduce(*d.histogram, comm, 0, true);

  if(params.verbosity >= 1) std::cout << "Done" << std::endl;

  return sol_sum;
}

} // namespace som
