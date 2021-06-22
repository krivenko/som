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

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>

#include <triqs/utility/signal_handler.hpp>

#include "kernels/all.hpp"
#include "objective_function.hpp"
#include "fit_quality.hpp"
#include "solution_worker.hpp"
#include "som_core.hpp"

namespace som {

using std::to_string;
using namespace mpi;
using namespace triqs::arrays;
using namespace triqs::gfs;
using triqs::statistics::histogram;

void fatal_error(std::string const& message) {
  TRIQS_RUNTIME_ERROR << "som_core: " << message;
}

void warning(std::string const& message) {
  std::cout << "WARNING: " << message << std::endl;
}

template <typename... GfOpts>
void check_input_gf(gf_const_view<GfOpts...> g, gf_const_view<GfOpts...> S) {
  if(g.mesh() != S.mesh() || g.target_shape() != S.target_shape())
    fatal_error(
        "input quantity and the error-bar function S must have equivalent "
        "structure");

  auto shape = g.target_shape();
  if(shape[0] != shape[1])
    fatal_error("matrix-valued input quantities must be square");
}

template <typename MeshType>
void check_gf_dim(gf_const_view<MeshType> g, int expected_dim) {
  auto shape = g.target_shape();
  if(shape[0] != expected_dim || shape[1] != expected_dim)
    fatal_error("expected a " + mesh_traits<MeshType>::name() +
                " Green's function with matrix dimensions " +
                to_string(expected_dim) + "x" + to_string(expected_dim));
}
template <typename MeshType>
void check_gf_dim(gf_view<MeshType> g, int expected_dim) {
  check_gf_dim(make_const_view(g), expected_dim);
}

template <typename MeshType>
void check_gf_stat(gf_const_view<MeshType> g, statistic_enum expected_stat) {
  if(g.domain().statistic != expected_stat)
    fatal_error("expected a " + mesh_traits<MeshType>::name() +
                " Green's function with " +
                (expected_stat == Fermion ? "fermionic" : "bosonic") +
                " statistics");
}
template <typename MeshType>
void check_gf_stat(gf_view<MeshType> g, statistic_enum expected_stat) {
  return check_gf_stat(make_const_view(g), expected_stat);
}

//////////////////////
// som_core::data_t //
//////////////////////

template <typename Mesh>
som_core::data_t::data_t(Mesh const& mesh, cache_index & ci) :
  rhs(input_data_t<Mesh>(mesh.size())),
  error_bars(input_data_t<Mesh>(mesh.size())),
  final_solution(ci)
{}

template <typename... GfOpts>
void som_core::data_t::init_input(int i,
                                  gf_const_view<GfOpts...> g,
                                  gf_const_view<GfOpts...> S,
                                  double norm_) {
  using mesh_t = typename gf_const_view<GfOpts...>::mesh_t;
  std::get<input_data_t<mesh_t>>(rhs)() = g.data()(range(), i, i);
  std::get<input_data_t<mesh_t>>(error_bars)() = S.data()(range(), i, i);
  norm = norm_;
}

////////////////////////////
// som_core: Constructors //
////////////////////////////

// Imaginary time
som_core::som_core(gf_const_view<imtime> g_tau, gf_const_view<imtime> S_tau,
                   observable_kind kind, vector<double> const& norms)
   : kind(kind)
   , mesh(g_tau.mesh())
   , data(g_tau.target_shape()[0], data_t(g_tau.mesh(), ci)) {

  if(is_stat_relevant(kind)) check_gf_stat(g_tau, observable_statistics(kind));

  check_input_gf(g_tau, S_tau);
  if(!is_gf_real(g_tau) || !is_gf_real(S_tau))
    fatal_error("imaginary time " + observable_name(kind) + " must be real");
  gf<imtime, matrix_real_valued> g_tau_real = real(g_tau),
                                 S_tau_real = real(S_tau);

  int gf_dim = g_tau.target_shape()[0];
  for(int i : range(gf_dim)) {
    data[i].init_input(i,
                       make_const_view(g_tau_real),
                       make_const_view(S_tau_real),
                       norms.size() > 0 ? norms[i] : 1.0);
  }
}

// Imaginary frequency
som_core::som_core(gf_const_view<imfreq> g_iw, gf_const_view<imfreq> S_iw,
                   observable_kind kind, vector<double> const& norms)
   : kind(kind)
   , mesh(g_iw.mesh())
   , data(g_iw.target_shape()[0], data_t(positive_freq_view(g_iw).mesh(), ci)) {

  if(is_stat_relevant(kind)) check_gf_stat(g_iw, observable_statistics(kind));

  check_input_gf(g_iw, S_iw);
  if(!is_gf_real_in_tau(g_iw) || !is_gf_real_in_tau(S_iw))
    fatal_error("imaginary frequency " + observable_name(kind) +
                R"( must be real in \tau-domain)");
  auto g_iw_pos = positive_freq_view(g_iw);
  auto S_iw_pos = positive_freq_view(S_iw);
  check_input_gf(g_iw_pos, S_iw_pos);

  int gf_dim = g_iw_pos.target_shape()[0];
  for(int i : range(gf_dim)) {
    data[i].init_input(i,
                       g_iw_pos,
                       S_iw_pos,
                       norms.size() > 0 ? norms[i] : 1.0);
  }
}

// Legendre coefficients
som_core::som_core(gf_const_view<legendre> g_l, gf_const_view<legendre> S_l,
                   observable_kind kind, vector<double> const& norms)
   : kind(kind)
   , mesh(g_l.mesh())
   , data(g_l.target_shape()[0], data_t(g_l.mesh(), ci)) {

  if(is_stat_relevant(kind)) check_gf_stat(g_l, observable_statistics(kind));

  check_input_gf(g_l, S_l);
  if(!is_gf_real(g_l) || !is_gf_real(S_l))
    fatal_error("Legendre " + observable_name(kind) + " must be real");
  gf<legendre, matrix_real_valued> g_l_real = real(g_l), S_l_real = real(S_l);

  int gf_dim = g_l_real.target_shape()[0];
  for(int i : range(gf_dim)) {
    data[i].init_input(i,
                       make_const_view(g_l_real),
                       make_const_view(S_l_real),
                       norms.size() > 0 ? norms[i] : 1.0);
  }
}

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

///////////////////////////////////////
// triqs_gf_view_assign_delegation() //
///////////////////////////////////////

void fill_data(gf_view<imtime> g_tau, int i, vector<double> const& data) {
  g_tau.data()(range(), i, i) = data;
}

void fill_data(gf_view<imfreq> g_iw, int i, vector<dcomplex> const& data) {
  auto g_positive_freq = positive_freq_view(g_iw);
  g_positive_freq.data()(range(), i, i) = data;
  g_iw = make_gf_from_real_gf(make_const_view(g_positive_freq));
}

void fill_data(gf_view<legendre> g_l, int i, vector<double> const& data) {
  g_l.data()(range(), i, i) = data;
}

template <typename MeshType>
void triqs_gf_view_assign_delegation(gf_view<MeshType> g,
                                     som_core const& cont) {
  auto gf_dim = cont.data.size();
  check_gf_dim(g, gf_dim);
  if(is_stat_relevant(cont.kind))
    check_gf_stat(g, observable_statistics(cont.kind));

  g() = 0;
#define FILL_DATA_CASE(r, d, ok)                                               \
  case ok: {                                                                   \
    kernel<ok, MeshType> kern(g.mesh());                                       \
    for(int i : range(gf_dim)) {                                               \
      fill_data(g, i, kern(cont.data[i].final_solution));                      \
    }                                                                          \
    return;                                                                    \
  }
  switch(cont.kind) {
    BOOST_PP_SEQ_FOR_EACH(FILL_DATA_CASE, _, ALL_OBSERVABLES)
    default: fatal_error("unknown observable kind " + to_string(cont.kind));
  }
#undef FILL_DATA_CASE
}

template <>
void triqs_gf_view_assign_delegation<refreq>(gf_view<refreq> g_w,
                                             som_core const& cont) {
  auto gf_dim = cont.data.size();
  check_gf_dim(g_w, gf_dim);
  for(int i : range(gf_dim)) {
    back_transform(cont.kind,
                   cont.data[i].final_solution,
                   const_cast<som_core&>(cont).ci,
                   slice_target_to_scalar(g_w, i, i),
                   const_cast<som_core&>(cont).comm
                   );
  }
}

template void triqs_gf_view_assign_delegation<imtime>(gf_view<imtime>,
                                                      som_core const&);
template void triqs_gf_view_assign_delegation<imfreq>(gf_view<imfreq>,
                                                      som_core const&);
template void triqs_gf_view_assign_delegation<legendre>(gf_view<legendre>,
                                                        som_core const&);

//////////////////////////////
// som_core::compute_tail() //
//////////////////////////////

array<dcomplex, 3> som_core::compute_tail(int max_order) const {
  auto gf_dim = data.size();
  array<dcomplex, 3> tail(max_order + 1, gf_dim, gf_dim);
  for(int i : range(gf_dim)) {
    tail(range(), i, i) =
        som::compute_tail(kind, data[i].final_solution,
                          const_cast<som_core*>(this)->ci,
                          const_cast<som_core*>(this)->comm,
                          max_order);
  }
  return tail;
}

/////////////////////////
// som_core: Accessors //
/////////////////////////

configuration const& som_core::get_solution(int i) const {
  if(i >= data.size())
    fatal_error("Matrix element index " + to_string(i) + " out of bounds");
  return data[i].final_solution;
}

std::vector<configuration> som_core::get_solutions() const {
  std::vector<configuration> conf;
  conf.reserve(data.size());
  for(auto const& d : data)
    conf.emplace_back(d.final_solution);
  return conf;
}

std::optional<histogram> const& som_core::get_histogram(int i) const {
  if(i >= data.size())
    fatal_error("Matrix element index " + to_string(i) + " out of bounds");
  return data[i].histogram;
}

std::optional<std::vector<histogram>> som_core::get_histograms() const {
  if(data.back().histogram) {
    std::vector<histogram> histograms;
    histograms.reserve(data.size());
    for(auto const& d : data)
      histograms.emplace_back(*d.histogram);
    return std::move(histograms);
  } else
    return {};
}

///////////////////////
// som_core::clear() //
///////////////////////

void som_core::clear() {
  for(auto & d : data) {
    d.basis_solutions.clear();
    d.final_solution.clear();
    d.histogram.reset();
  }
}

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
                      << ", D = " << worker.get_objf_value() << ")." << std::endl;
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
    fatal_error("Wrong f_range in adjust_f()");

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
