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

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/seq/enum.hpp>

#include <mpi/mpi.hpp>
#include <nda/nda.hpp>

#include <triqs/utility/signal_handler.hpp>

#include <som/kernels/all.hpp>
#include <som/solution_worker.hpp>

#include "common.hxx"
#include "som_core.hpp"

namespace som {

using triqs::stat::histogram;

void accumulate_parameters_t::validate(observable_kind kind) const {
  worker_parameters_t::validate(kind);

  using std::to_string;

  if(f <= 0)
    fatal_error("Number of global updates f must be positive (got f = " +
                to_string(f) + ")");

  if(l <= 0)
    fatal_error("Number of particular solutions l must be positive (got l = " +
                to_string(l) + ")");

  if(adjust_l) {
    if(adjust_l_range.first <= 0 || adjust_l_range.second <= 0 ||
       adjust_l_range.first > adjust_l_range.second)
      fatal_error("Wrong adjust_l_range = [" + to_string(adjust_l_range.first) +
                  ";" + to_string(adjust_l_range.second) + "]");

    if(adjust_l_good_chi < 1.0)
      fatal_error("Parameter adjust_l_good_chi must be >= 1.0 (got " +
                  to_string(adjust_l_good_chi) + ")");

    if(adjust_l_verygood_chi < 1.0)
      fatal_error("Parameter adjust_l_verygood_chi must be >= 1.0 (got " +
                  to_string(adjust_l_verygood_chi) + ")");

    if(adjust_l_verygood_chi > adjust_l_good_chi)
      fatal_error("adjust_l_verygood_chi cannot exceed adjust_l_good_chi");

    if(adjust_l_ratio < 0 || adjust_l_ratio > 1.0)
      fatal_error("Parameter adjust_l_ratio = " + to_string(adjust_l_ratio) +
                  " is not in [0;1.0]");
  }

  if(hist_max <= 1.0)
    fatal_error("Parameter hist_max = " + to_string(hist_max) +
                " must exceed 1.0");

  if(make_histograms && hist_n_bins < 1)
    fatal_error("Histograms must contain at least one bin (got hist_n_bins = " +
                to_string(hist_n_bins) + ")");
}

////////////////////////////
// som_core::accumulate() //
////////////////////////////

template <typename KernelType> void som_core::accumulate_impl() {

  auto stop_callback = triqs::utility::clock_callback(params.max_time);

  auto const& m = std::get<typename KernelType::mesh_type>(mesh);

  if(params.verbosity > 0) {
    mpi_cout(comm) << "Constructing integral kernel... " << std::endl;
  }
  KernelType kernel(m);
  if(params.verbosity > 0) {
    mpi_cout(comm) << "Constructed kernel: " << kernel << std::endl;
  }

  // Find solution for each component of GF
  for(int n = 0; n < data.size(); ++n) {
    auto& d = data[n];

    if(params.verbosity > 0)
      mpi_cout(comm)
          << "Accumulating particular solutions for observable component [" << n
          << "," << n << "]" << std::endl;

    auto of = d.make_objf(kernel);
    if(params.verbosity >= 2 && of.get_U_dagger()) {
      mpi_cout(comm) << "Supplied covariance matrix has " +
                            std::to_string(of.get_sigma2().size())
                     << " positive eigenvalues" << std::endl;
    }

    solution_worker<KernelType> worker(
        of, d.norm, ci, params, stop_callback, params.f);
    auto& rng = worker.get_rng();

    // Reset final solution as it is no more valid
    d.final_solution.clear();

    int n_sol_max = 0;    // Number of solutions to be accumulated
    int n_sol = 0, i = 0; // Global and rank-local indices of solution
    int n_good_solutions = 0,
        n_verygood_solutions = 0; // Number of good and very good solutions
    do {
      if(params.adjust_l) {
        n_sol_max += params.adjust_l_range.first;
        if(n_sol_max > params.adjust_l_range.second) {
          if(params.verbosity >= 1)
            mpi_cout(comm)
                << "WARNING: Upper bound of adjust_l_range has been reached"
                << std::endl;
          break;
        }
      } else
        n_sol_max = params.l;

      d.particular_solutions.reserve(d.particular_solutions.size() + n_sol_max);

      if(params.verbosity >= 1) {
        if(params.adjust_l) {
          mpi_cout(comm) << "Increasing the total number of solutions to be "
                            "accumulated to "
                         << n_sol_max << std::endl;
        } else {
          mpi_cout(comm) << "Accumulating a total of " << n_sol_max
                         << " solutions" << std::endl;
        }
      }

      for(; (n_sol = comm.rank() + i * comm.size()) < n_sol_max; ++i) {
        if(params.verbosity >= 2) {
          mpi_cout(comm) << "Accumulation of particular solution " << n_sol
                         << std::endl;
        }

        d.particular_solutions.emplace_back(worker(1 + rng(params.max_rects)),
                                            0);

        double chi2 = worker.get_objf_value();
        d.particular_solutions.back().second = chi2;
        d.objf_min = std::min(d.objf_min, chi2);

        if(params.verbosity >= 2) {
          mpi_cout(comm) << "Solution " << n_sol << ": χ = " << std::sqrt(chi2)
                         << std::endl;
        }
      }
      comm.barrier(0);

      // Global minimum of \chi^2_min
      // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
      d.objf_min = mpi::all_reduce(d.objf_min, comm, MPI_MIN);
      double chi_min = std::sqrt(d.objf_min);

      // Recalculate numbers of good and very good solutions
      n_good_solutions = n_verygood_solutions = 0;
      for(auto const& s : d.particular_solutions) {
        if(std::sqrt(s.second) / chi_min <= params.adjust_l_good_chi)
          ++n_good_solutions;
        if(std::sqrt(s.second) / chi_min <= params.adjust_l_verygood_chi)
          ++n_verygood_solutions;
      }
      n_good_solutions = mpi::all_reduce(n_good_solutions, comm);
      n_verygood_solutions = mpi::all_reduce(n_verygood_solutions, comm);

      if(params.verbosity >= 1) {
        mpi_cout(comm) << "χ_{min} = " << chi_min << std::endl;
        if(params.adjust_l) {
          mpi_cout(comm) << "Number of good solutions (χ / χ_{min} ≤ "
                         << params.adjust_l_good_chi
                         << ") = " << n_good_solutions << std::endl;
          mpi_cout(comm) << "Number of very good solutions (χ / χ_{min} ≤ "
                         << params.adjust_l_verygood_chi
                         << ") = " << n_verygood_solutions << std::endl;
        }
      }

    } while(params.adjust_l &&
            double(n_verygood_solutions) / double(n_good_solutions) <
                params.adjust_l_ratio);

    comm.barrier(0);

    // Recompute the histograms
    if(params.make_histograms) {
      d.histogram = histogram(std::sqrt(d.objf_min),
                              std::sqrt(d.objf_min) * params.hist_max,
                              params.hist_n_bins);
      histogram& h = *d.histogram;
      for(auto const& s : d.particular_solutions) h << std::sqrt(s.second);
      h = mpi::all_reduce(h, comm);
    }

    if(params.verbosity >= 1) {
      mpi_cout(comm) << "Accumulation completed" << std::endl;
    }
  }

  ci.invalidate_all();

  if(params.verbosity >= 1)
    mpi_cout(comm)
        << "Accumulation for all observable components has been completed"
        << std::endl;
}

void som_core::accumulate(accumulate_parameters_t const& p) {
  params = p;
  params.validate(kind);

  triqs::signal_handler::start();
  accumulate_status = 0;
  try {
#define IMPL_CASE(r, okmk)                                                     \
  case(kernel_id<BOOST_PP_SEQ_ELEM(1, okmk)>(BOOST_PP_SEQ_ELEM(0, okmk))):     \
    accumulate_impl<kernel<BOOST_PP_SEQ_ENUM(okmk)>>();                        \
    break;

    SELECT_KERNEL(IMPL_CASE, som_core::accumulate())
#undef IMPL_CASE
  } catch(stopped& e) {
    accumulate_status = e.code;
    triqs::signal_handler::received(true);
  }
  triqs::signal_handler::stop();
}

} // namespace som
