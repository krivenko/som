/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2025 Igor Krivenko
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
#include <string>

#include <boost/preprocessor/seq/for_each_product.hpp>

#include <triqs/utility/signal_handler.hpp>

#include <som/kernels/all.hpp>
#include <som/solution_functionals/fit_quality.hpp>
#include <som/solution_functionals/objective_function.hpp>
#include <som/solution_worker.hpp>

#include "common.hxx"
#include "som_core.hpp"

namespace som {

using namespace triqs::gfs;

void adjust_f_parameters_t::validate(observable_kind kind) const {
  worker_parameters_t::validate(kind);

  using std::to_string;

  if(f_range.first <= 0 || f_range.second <= 0 ||
     f_range.first > f_range.second)
    fatal_error("Wrong f_range = [" + to_string(f_range.first) + ";" +
                to_string(f_range.second) + "]");

  if(l <= 0)
    fatal_error("Number of particular solutions l must be positive (got l = " +
                to_string(l) + ")");

  if(kappa <= 0 || kappa > 0.5)
    fatal_error("Parameter kappa = " + to_string(kappa) + " is not in ]0;0.5]");
}

//////////////////////////
// som_core::adjust_f() //
//////////////////////////

template <typename KernelType>
int som_core::adjust_f_impl(adjust_f_parameters_t const& p) {

  int F_max = p.f_range.first;

  auto const& m = std::get<typename KernelType::mesh_type>(mesh);

  triqs::signal_handler::start();
  try {
    auto stop_callback = triqs::utility::clock_callback(p.max_time);

    if(p.verbosity > 0) {
      mpi_cout(comm) << "Constructing integral kernel... \n";
    }
    KernelType kernel(m);
    if(p.verbosity > 0) {
      mpi_cout(comm) << "Constructed kernel: " << kernel << '\n';
    }

    // Find solution for each component of GF
    for(int n = 0; n < data.size(); ++n) {
      auto& d = data[n];

      if(p.verbosity > 0)
        mpi_cout(comm) << "Running algorithm for observable component [" << n
                       << "," << n << "]\n";

      auto of = d.make_objf<KernelType>(kernel);
      if(of.get_U_dagger())
        fatal_error(
            "Automatic adjustment procedure for the number of global updates "
            "(F) is not available when a full covariance matrix of the RHS is "
            "used");

      typename fit_quality<KernelType>::rhs_type error_bars =
          sqrt(of.get_sigma2());
      fit_quality<KernelType> fq(kernel, of.get_rhs(), error_bars);

      int F = F_max;

      int l_good = 0;
      for(;; l_good = 0, F *= 2) {
        // Upper bound of f_range is reached
        if(F >= p.f_range.second) {
          F = p.f_range.second;
          if(p.verbosity >= 1)
            mpi_cout(comm)
                << "WARNING: Upper bound of f_range has been reached,"
                   " will use F = " +
                       std::to_string(F)
                << '\n';
          break;
        }

        solution_worker<KernelType> worker(of, d.norm, ci, p, stop_callback, F);
        // cppcheck-suppress constVariableReference
        auto& rng = worker.get_rng();

        int n_sol = {};
        for(int i = 0; (n_sol = comm.rank() + i * comm.size()) < p.l; ++i) {
          if(p.verbosity >= 2) {
            mpi_cout(comm) << "Accumulation of particular solution " << n_sol
                           << '\n';
          }

          auto solution = worker(1 + rng(params.max_rects));
          double kappa = fq(solution);

          if(kappa > p.kappa) ++l_good;
          if(p.verbosity >= 2) {
            mpi_cout(comm) << "Particular solution " << n_sol << " is "
                           << (kappa > p.kappa ? "" : "not ") << R"(good (κ = )"
                           << kappa
                           << ", χ = " << std::sqrt(worker.get_objf_value())
                           << ").\n";
          }
        }
        comm.barrier(0);
        l_good = mpi::all_reduce(l_good, comm);

        if(p.verbosity >= 1)
          mpi_cout(comm) << "F = " << F << ", " << l_good
                         << R"( solutions with κ > )" << p.kappa << " (out of "
                         << p.l << ")\n";

        // Converged
        if(l_good > p.l / 2) {
          if(p.verbosity >= 1) mpi_cout(comm) << "F = " << F << " is enough.\n";
          break;
        }
      }

      F_max = std::max(F_max, F);
    }
  } catch(stopped& e) { triqs::signal_handler::received(true); }

  ci.invalidate_all();

  triqs::signal_handler::stop();

  return F_max;
}

int som_core::adjust_f(adjust_f_parameters_t const& p) {
  p.validate(kind);

  if(p.verbosity >= 1) {
    mpi_cout(comm) << "Adjusting the number of global updates F using " << p.l
                   << " particular solutions ...\n";
  }

#define IMPL_CASE(r, okmk)                                                     \
  case(kernel_id<BOOST_PP_SEQ_ELEM(1, okmk)>(BOOST_PP_SEQ_ELEM(0, okmk))):     \
    return adjust_f_impl<kernel<BOOST_PP_SEQ_ENUM(okmk)>>(p);

  SELECT_KERNEL(IMPL_CASE, som_core::adjust_f())
#undef IMPL_CASE
}

} // namespace som
