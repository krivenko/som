/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2018 by I. Krivenko
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
#include <limits>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>

#include "som_core.hpp"
#include "kernels/fermiongf_imtime.hpp"
#include "kernels/fermiongf_imfreq.hpp"
#include "kernels/fermiongf_legendre.hpp"
#include "kernels/bosoncorr_imtime.hpp"
#include "kernels/bosoncorr_imfreq.hpp"
#include "kernels/bosoncorr_legendre.hpp"
#include "kernels/bosonautocorr_imtime.hpp"
#include "kernels/bosonautocorr_imfreq.hpp"
#include "kernels/bosonautocorr_legendre.hpp"
#include "kernels/zerotemp_imtime.hpp"
#include "kernels/zerotemp_imfreq.hpp"
#include "kernels/zerotemp_legendre.hpp"
#include "objective_function.hpp"
#include "fit_quality.hpp"
#include "solution_worker.hpp"

namespace som {

using namespace triqs::mpi;
using std::to_string;

void fatal_error(std::string const& message) {
 TRIQS_RUNTIME_ERROR << "som_core: " << message;
}

void warning(std::string const& message) {
 std::cout << "WARNING: " << message << std::endl;
}

template<typename... GfOpts>
void check_input_gf(gf_const_view<GfOpts...> g, gf_const_view<GfOpts...> S) {
 if(g.mesh() != S.mesh() || get_target_shape(g) != get_target_shape(S))
  fatal_error("input quantity and the error-bar function S must have equivalent structure");

 auto shape = get_target_shape(g);
 if(shape[0] != shape[1]) fatal_error("matrix-valued input quantities must be square");
}

template<typename... GfOpts>
void som_core::set_input_data(gf_const_view<GfOpts...> g, gf_const_view<GfOpts...> S) {
 using mesh_t = typename gf_const_view<GfOpts...>::mesh_t;

 auto & rhs_ = (input_data_t<mesh_t>&)rhs;
 auto & error_bars_ = (input_data_t<mesh_t>&)error_bars;
 auto gf_dim = get_target_shape(g)[0];

 results.reserve(gf_dim);
 for(int i = 0; i < gf_dim; ++i) {
  rhs_.emplace_back(g.data()(range(),i,i));
  error_bars_.emplace_back(S.data()(range(),i,i));
  results.emplace_back(ci);
 }
 histograms.resize(gf_dim);
}

template<typename MeshType> void check_gf_dim(gf_const_view<MeshType> g, int expected_dim) {
 auto shape = get_target_shape(g);
 if(shape[0] != expected_dim || shape[1] != expected_dim)
  fatal_error("expected a " + mesh_traits<MeshType>::name()
              + " Green's function with matrix dimensions "
              + to_string(expected_dim) + "x" + to_string(expected_dim));
}
template<typename MeshType> void check_gf_dim(gf_view<MeshType> g, int expected_dim) {
 check_gf_dim(make_const_view(g), expected_dim);
}

template<typename MeshType> void check_gf_stat(gf_const_view<MeshType> g,
                                               triqs::gfs::statistic_enum expected_stat) {
 if(g.domain().statistic != expected_stat)
  fatal_error("expected a " + mesh_traits<MeshType>::name()
              + " Green's function with "
              + (expected_stat == Fermion ? "fermionic" : "bosonic")
              + " statistics");
}
template<typename MeshType> void check_gf_stat(gf_view<MeshType> g,
                                               triqs::gfs::statistic_enum expected_stat) {
 return check_gf_stat(make_const_view(g), expected_stat);
}

vector<double> make_default_norms(vector<double> const& norms, int dim) {
 if(norms.size() != 0) return norms;
 else {
  vector<double> def_norms(dim);
  def_norms() = 1.0;
  return def_norms;
 }
}

//////////////////
// Constructors //
//////////////////

// Imaginary time
som_core::som_core(gf_const_view<imtime> g_tau, gf_const_view<imtime> S_tau,
                   observable_kind kind, vector<double> const& norms) :
 mesh(g_tau.mesh()), kind(kind), norms(make_default_norms(norms,get_target_shape(g_tau)[0])),
 rhs(input_data_r_t()), error_bars(input_data_r_t()) {

 if(is_stat_relevant(kind)) check_gf_stat(g_tau, observable_statistics(kind));

 check_input_gf(g_tau,S_tau);
 if(!is_gf_real(g_tau) || !is_gf_real(S_tau))
  fatal_error("imaginary time " + observable_name(kind) + " must be real");
 gf<imtime, matrix_real_valued> g_tau_real = real(g_tau), S_tau_real = real(S_tau);
 set_input_data(make_const_view(g_tau_real), make_const_view(S_tau_real));
}

// Imaginary frequency
som_core::som_core(gf_const_view<imfreq> g_iw, gf_const_view<imfreq> S_iw,
                   observable_kind kind, vector<double> const& norms) :
 mesh(g_iw.mesh()), kind(kind), norms(make_default_norms(norms,get_target_shape(g_iw)[0])),
 rhs(input_data_c_t()), error_bars(input_data_c_t()) {

 if(is_stat_relevant(kind)) check_gf_stat(g_iw, observable_statistics(kind));

 check_input_gf(g_iw,S_iw);
 if(!is_gf_real_in_tau(g_iw) || !is_gf_real_in_tau(S_iw))
  fatal_error("imaginary frequency " + observable_name(kind) + " must be real in \\tau-domain");
 auto g_iw_pos = positive_freq_view(g_iw);
 auto S_iw_pos = positive_freq_view(S_iw);
 check_input_gf(g_iw_pos,S_iw_pos);
 set_input_data(g_iw_pos,S_iw_pos);
}

// Legendre coefficients
som_core::som_core(gf_const_view<legendre> g_l, gf_const_view<legendre> S_l,
                   observable_kind kind, vector<double> const& norms) :
 mesh(g_l.mesh()), kind(kind), norms(make_default_norms(norms,get_target_shape(g_l)[0])),
 rhs(input_data_r_t()), error_bars(input_data_r_t()) {

 if(is_stat_relevant(kind)) check_gf_stat(g_l, observable_statistics(kind));

 check_input_gf(g_l,S_l);
 if(!is_gf_real(g_l) || !is_gf_real(S_l))
  fatal_error("Legendre " + observable_name(kind) + " must be real");
 gf<legendre, matrix_real_valued> g_l_real = real(g_l), S_l_real = real(S_l);
 set_input_data(make_const_view(g_l_real), make_const_view(S_l_real));
}

///////////
// run() //
///////////

void som_core::run(run_parameters_t const& p) {

 params = p;

 double e_min, e_max;
 std::tie(e_min,e_max) = max_energy_window(kind);
 if(params.energy_window.first < e_min) {
  params.energy_window.first = e_min;
  if(params.verbosity > 0)
   warning("left boundary of the energy window is reset to " + std::to_string(e_min));
 }
 if(params.energy_window.second > e_max) {
  params.energy_window.second = e_max;
  if(params.verbosity > 0)
   warning("right boundary of the energy window is reset to " + std::to_string(e_max));
 }

 if(params.energy_window.first >= params.energy_window.second)
  fatal_error("wrong energy window [" + to_string(params.energy_window.first) +
              ";" + to_string(params.energy_window.second) + "]");

 double window_width = params.energy_window.second - params.energy_window.first;
 if(params.min_rect_width <= 0 || params.min_rect_width > 1)
  fatal_error("min_rect_width must be in (0;1]");

 if(params.min_rect_weight <= 0 || params.min_rect_weight > 1)
  fatal_error("min_rect_weight must be in (0;1]");

 if(params.adjust_f_range.first > params.adjust_f_range.second)
  fatal_error("Wrong adjust_f_range");

 if(params.adjust_l_range.first > params.adjust_l_range.second)
  fatal_error("Wrong adjust_l_range");

 triqs::signal_handler::start();
 run_status = 0;
 try {
  #define RUN_IMPL_CASE(r, okmk)                                              \
   case (int(BOOST_PP_SEQ_ELEM(0,okmk)) +                                     \
         n_observable_kinds * mesh_traits<BOOST_PP_SEQ_ELEM(1,okmk)>::index): \
         run_impl<kernel<BOOST_PP_SEQ_ENUM(okmk)>>(); break;
  switch(int(kind) + n_observable_kinds * mesh.index()) {
   BOOST_PP_SEQ_FOR_EACH_PRODUCT(RUN_IMPL_CASE, (ALL_OBSERVABLES)(ALL_INPUT_MESHES))
  }
  #undef RUN_IMPL_CASE
 } catch(stopped & e) {
  run_status = e.code;
  triqs::signal_handler::received(true);
 }
 triqs::signal_handler::stop();
}

template<typename KernelType> void som_core::run_impl() {

 auto stop_callback = triqs::utility::clock_callback(params.max_time);

 using mesh_t = typename KernelType::mesh_type;
 mesh_t const& m = mesh;

 if(params.verbosity > 0) {
  std::cout << "Constructing integral kernel... " << std::flush;
 }
 KernelType kernel(m);
 if(params.verbosity > 0) {
  std::cout << "done" << std::endl << std::flush;
  std::cout << "Kernel: " << kernel << std::endl;
 }

 // Find solution for each component of GF
 for(int i = 0; i < results.size(); ++i) {
  if(params.verbosity > 0)
   std::cout << "Running algorithm for observable component [" << i << "," << i << "]" << std::endl;

  auto const& rhs_ = (input_data_t<mesh_t> const&) rhs;
  auto const& error_bars_ = (input_data_t<mesh_t> const&) error_bars;
  double norm = norms[i];

  int F = params.adjust_f ? adjust_f(kernel,rhs_[i],error_bars_[i],norm,stop_callback) : params.f;
  results[i] = accumulate(kernel,rhs_[i],error_bars_[i],norm,histograms[i],stop_callback,F);
 }

 ci.invalidate_all();
}

template<typename KernelType> int som_core::adjust_f(KernelType const& kern,
 typename KernelType::result_type rhs_,
 typename KernelType::result_type error_bars_,
 double norm,
 std::function<bool()> const& stop_callback) {

 if(params.verbosity >= 1) {
  std::cout << "Adjusting the number of global updates F using "
            << params.adjust_f_l << " particular solutions ..." << std::endl;
 }
 objective_function<KernelType> of(kern, rhs_, error_bars_);
 fit_quality<KernelType> fq(kern, rhs_, error_bars_);

 int F = params.adjust_f_range.first;

 int l_good;
 for(l_good = 0;;l_good = 0, F *= 2) {
  // Upper bound of adjust_f_range is reached
  if(F >= params.adjust_f_range.second) {
   F = params.adjust_f_range.second;
   if(params.verbosity >= 1)
    warning("Upper bound of adjust_f_range has been reached, will use F = " + std::to_string(F));
   break;
  }

  solution_worker<KernelType> worker(of,norm,ci,params,stop_callback,F);
  auto & rng = worker.get_rng();

  int n_sol;
  for(int i = 0; (n_sol = comm.rank() + i*comm.size()) < params.adjust_f_l; ++i) {
   if(params.verbosity >= 2) {
    std::cout << "[Node " << comm.rank() << "] Accumulation of particular solution "
              << n_sol << std::endl;
   }

   auto solution = worker(1 + rng(params.max_rects));
   double kappa = fq(solution);

   if(kappa > params.adjust_f_kappa) ++l_good;
   if(params.verbosity >= 2) {
    std::cout << "[Node " << comm.rank() << "] Particular solution "
              << n_sol << " is " << (kappa > params.adjust_f_kappa ? "" : "not ")
              << "good (\\kappa = " << kappa
              << ", D = " << worker.get_objf_value() << ")." << std::endl;
   }
  }
  comm.barrier();
  l_good = mpi_all_reduce(l_good,comm);

  if(params.verbosity >= 1)
   std::cout << "F = " << F << ", "
             << l_good<< " solutions with \\kappa > " << params.adjust_f_kappa
             << " (out of " << params.adjust_f_l << ")" << std::endl;

  // Converged
  if(l_good > params.adjust_f_l/2) {
   if(params.verbosity >= 1) std::cout << "F = " << F << " is enough." << std::endl;
   break;
  }
 }

 return F;
}

template<typename KernelType> configuration som_core::accumulate(KernelType const& kern,
 typename KernelType::result_type rhs_,
 typename KernelType::result_type error_bars_,
 double norm,
 histogram & hist,
 std::function<bool()> const& stop_callback,
 int F) {

 if(params.verbosity >= 1) std::cout << "Accumulating particular solutions ..." << std::endl;

 objective_function<KernelType> of(kern, rhs_, error_bars_);
 solution_worker<KernelType> worker(of,norm,ci,params,stop_callback,F);
 auto & rng = worker.get_rng();

 // Pairs (configuration, objective function)
 std::vector<std::pair<configuration,double>> solutions;

 int n_sol_max = 0;                          // Number of solutions to be accumulated
 int n_sol, i = 0;                           // Global and rank-local indices of solution
 int n_good_solutions, n_verygood_solutions; // Number of good and very good solutions
 double objf_min = HUGE_VAL;                 // Minimal value of D
 do {
  if(params.adjust_l) {
   n_sol_max += params.adjust_l_range.first;
   if(n_sol_max > params.adjust_l_range.second) {
    if(params.verbosity >= 1) warning("Upper bound of adjust_l_range has been reached");
    break;
   }
  } else
   n_sol_max = params.l;
  solutions.reserve(n_sol_max);

  if(params.verbosity >= 1)
   std::cout << "Increasing the total number of solutions to be accumulated to "
             << n_sol_max << std::endl;

  for(; (n_sol = comm.rank() + i*comm.size()) < n_sol_max; ++i) {
   if(params.verbosity >= 2) {
    std::cout << "[Node " << comm.rank() << "] Accumulation of particular solution "
              << n_sol << std::endl;
   }

   solutions.emplace_back(worker(1 + rng(params.max_rects)), 0);

   double D = worker.get_objf_value();
   solutions.back().second = D;
   objf_min = std::min(objf_min, D);

   if(params.verbosity >= 2) {
    std::cout << "[Node " << comm.rank() << "] Solution " << n_sol
              << ": D = " << D << std::endl;
   }
  }
  comm.barrier();

  // Global minimum of D_min
  objf_min = mpi_all_reduce(objf_min, comm, 0, MPI_MIN);

  // Recalculate numbers of good and very good solutions
  n_good_solutions = n_verygood_solutions = 0;
  for(auto const& s : solutions) {
   if(s.second/objf_min <= params.adjust_l_good_d) ++n_good_solutions;
   if(s.second/objf_min <= params.adjust_l_verygood_d) ++n_verygood_solutions;
  }
  n_good_solutions = mpi_all_reduce(n_good_solutions);
  n_verygood_solutions = mpi_all_reduce(n_verygood_solutions);

  if(params.verbosity >= 1) {
   std::cout << "D_min = " << objf_min << std::endl;
   std::cout << "Number of good solutions (D/D_min <= "
             << params.adjust_l_good_d << ") = " << n_good_solutions << std::endl;
   std::cout << "Number of very good solutions (D/D_min <= "
             << params.adjust_l_verygood_d << ") = " << n_verygood_solutions << std::endl;
  }

 } while(params.adjust_l &&
         double(n_verygood_solutions) / double(n_good_solutions) < params.adjust_l_ratio);

 comm.barrier();

 if(params.verbosity >= 1) {
  std::cout << "Accumulation complete." << std::endl;
  std::cout << "Summing up good solutions ..." << std::endl;
 }

 if(params.make_histograms)
  hist = histogram(objf_min, objf_min * params.hist_max, params.hist_n_bins);

 configuration sol_sum(ci);

 // Rank-local stage of summation
 for(auto const& s : solutions) {
  if(params.make_histograms) hist << s.second;
  // Pick only good solutions
  if(s.second/objf_min <= params.adjust_l_good_d) sol_sum += s.first;
 }
 sol_sum *= 1.0/double(n_good_solutions);

 // Sum over all processes
 sol_sum = mpi_reduce(sol_sum, comm, 0, true);

 if(params.make_histograms) hist = mpi_reduce(hist, comm, 0, true);

 if(params.verbosity >= 1) std::cout << "Done" << std::endl;

 return sol_sum;
}

///////////////////////////////////////
// triqs_gf_view_assign_delegation() //
///////////////////////////////////////

void fill_data(gf_view<imtime> g_tau, int i, vector<double> const& data) {
 g_tau.data()(range(),i,i) = data;
}

void fill_data(gf_view<imfreq> g_iw, int i, vector<dcomplex> const& data) {
 auto g_positive_freq = positive_freq_view(g_iw);
 g_positive_freq.data()(range(),i,i) = data;
 g_iw = make_gf_from_real_gf(make_const_view(g_positive_freq));
}

void fill_data(gf_view<legendre> g_l, int i, vector<double> const& data) {
 g_l.data()(range(),i,i) = data;
}

template<typename MeshType>
void triqs_gf_view_assign_delegation(gf_view<MeshType> g, som_core const& cont) {
 auto gf_dim = cont.results.size();
 check_gf_dim(g, gf_dim);
 if(is_stat_relevant(cont.kind)) check_gf_stat(g, observable_statistics(cont.kind));

 g() = 0;
#define FILL_DATA_CASE(r, data, ok)                                  \
 case ok: {                                                          \
  kernel<ok,MeshType> kern(g.mesh());                                \
  for(int i : range(gf_dim)) fill_data(g, i, kern(cont.results[i])); \
  return;                                                            \
 }
 switch(cont.kind) {
  BOOST_PP_SEQ_FOR_EACH(FILL_DATA_CASE, _, ALL_OBSERVABLES)
  default: fatal_error("unknown observable kind " + to_string(cont.kind));
 }
#undef FILL_DATA_CASE
}

template<> void triqs_gf_view_assign_delegation<refreq>(gf_view<refreq> g_w, som_core const& cont) {
 auto gf_dim = cont.results.size();
 check_gf_dim(g_w, gf_dim);
 for(int i : range(gf_dim))
  back_transform(cont.kind,
                 cont.results[i],
                 const_cast<som_core&>(cont).ci,
                 slice_target_to_scalar(g_w, i, i));
}

template void triqs_gf_view_assign_delegation<imtime>(gf_view<imtime>, som_core const&);
template void triqs_gf_view_assign_delegation<imfreq>(gf_view<imfreq>, som_core const&);
template void triqs_gf_view_assign_delegation<legendre>(gf_view<legendre>, som_core const&);

}
