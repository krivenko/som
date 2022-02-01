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

#include <cmath>
#include <iostream>

#include <boost/preprocessor/seq/for_each_product.hpp>

#include <triqs/utility/signal_handler.hpp>

#include <som/kernels/all.hpp>
#include <som/solution_functionals/objective_function.hpp>
#include <som/solution_functionals/fit_quality.hpp>
#include <som/solution_worker.hpp>

#include "som_core.hpp"

namespace som {

using namespace triqs::gfs;

//////////////////////////
// som_core::adjust_f() //
//////////////////////////

template <typename KernelType>
int som_core::adjust_f_impl(adjust_f_parameters_t const& p) {

  int F_max = p.f_range.first;

  using mesh_t = typename KernelType::mesh_type;
  mesh_t const& m = std::get<mesh_t>(mesh);

  triqs::signal_handler::start();
  try {
    auto stop_callback = triqs::utility::clock_callback(p.max_time);

    if(p.verbosity > 0) {
      std::cout << "Constructing integral kernel... " << std::flush;
    }
    KernelType kernel(m);
    if(p.verbosity > 0) {
      std::cout << "done" << std::endl;
      std::cout << "Kernel: " << kernel << std::endl;
    }

    // Find solution for each component of GF
    for(int n = 0; n < data.size(); ++n) {
      auto & d = data[n];

      if(p.verbosity > 0)
        std::cout << "Running algorithm for observable component [" << n << ","
                  << n << "]" << std::endl;

      auto const& rhs = d.get_rhs<mesh_t>();
      auto const& error_bars = d.get_error_bars<mesh_t>();

      objective_function<KernelType> of(kernel, rhs, error_bars);
      fit_quality<KernelType> fq(kernel, rhs, error_bars);

      int F = F_max;

      int l_good;
      for(l_good = 0;; l_good = 0, F *= 2) {
        // Upper bound of f_range is reached
        if(F >= p.f_range.second) {
          F = p.f_range.second;
          if(p.verbosity >= 1)
            std::cout << "WARNING: Upper bound of f_range has been reached,"
                         " will use F = " + std::to_string(F) << std::endl;
          break;
        }

        solution_worker<KernelType> worker(of,
                                          d.norm,
                                          ci,
                                          p,
                                          stop_callback,
                                          F);
        auto& rng = worker.get_rng();

        int n_sol;
        for(int i = 0; (n_sol = comm.rank() + i * comm.size()) < p.l;
            ++i) {
          if(p.verbosity >= 2) {
            std::cout << "[Rank " << comm.rank()
                      << "] Accumulation of particular solution " << n_sol
                      << std::endl;
          }

          auto solution = worker(1 + rng(params.max_rects));
          double kappa = fq(solution);

          if(kappa > p.kappa) ++l_good;
          if(p.verbosity >= 2) {
            std::cout << "[Rank " << comm.rank() << "] Particular solution "
                      << n_sol << " is "
                      << (kappa > p.kappa ? "" : "not ")
                      << R"(good (\kappa = )" << kappa
                      << ", \\chi = "
                      << std::sqrt(worker.get_objf_value()) << ")." << std::endl;
          }
        }
        comm.barrier();
        l_good = mpi::all_reduce(l_good, comm);

        if(p.verbosity >= 1)
          std::cout << "F = " << F << ", " << l_good
                    << R"( solutions with \kappa > )" << p.kappa
                    << " (out of " << p.l << ")" << std::endl;

        // Converged
        if(l_good > p.l / 2) {
          if(p.verbosity >= 1)
            std::cout << "F = " << F << " is enough." << std::endl;
          break;
        }
      }

      F_max = std::max(F_max, F);
    }
  } catch(stopped& e) {
    triqs::signal_handler::received(true);
  }

  ci.invalidate_all();

  triqs::signal_handler::stop();

  return F_max;
}

int som_core::adjust_f(adjust_f_parameters_t const& p) {
  if(p.f_range.first > p.f_range.second)
    TRIQS_RUNTIME_ERROR << "som_core: Wrong f_range in adjust_f()";

  if(p.verbosity >= 1) {
    std::cout << "Adjusting the number of global updates F using "
              << p.l << " particular solutions ..." << std::endl;
  }

#define RUN_IMPL_CASE(r, okmk)                                                 \
  case(int(BOOST_PP_SEQ_ELEM(0, okmk)) +                                       \
       n_observable_kinds * mesh_traits<BOOST_PP_SEQ_ELEM(1, okmk)>::index):   \
    return adjust_f_impl<kernel<BOOST_PP_SEQ_ENUM(okmk)>>(p);
  switch(int(kind) + n_observable_kinds * mesh.index()) {
    BOOST_PP_SEQ_FOR_EACH_PRODUCT(RUN_IMPL_CASE,
                                  (ALL_OBSERVABLES)(ALL_INPUT_MESHES))
    default: TRIQS_RUNTIME_ERROR << "som_core: unknown observable kind "
                                 << std::to_string(kind)
                                 << " in adjust_f()";
  }
#undef RUN_IMPL_CASE
}

} // namespace som
