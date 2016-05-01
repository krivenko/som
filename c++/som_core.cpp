/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016 by I. Krivenko
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include <limits>

#include "som_core.hpp"
#include "solution_worker.hpp"

namespace som {

using namespace triqs::mpi;
using std::to_string;

void fatal_error(std::string const& message) {
 TRIQS_RUNTIME_ERROR << "som_core: " << message;
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

//////////////////
// Constructors //
//////////////////

// Imaginary time
som_core::som_core(gf_const_view<imtime> g_tau, gf_const_view<imtime> S_tau,
                   observable_kind kind, double norm) :
 mesh(g_tau.mesh()), kind(kind), norm(norm),
 rhs(input_data_r_t()), error_bars(input_data_r_t()) {
 switch(kind) {
  case FermionGf: {
   if(g_tau.domain().statistic != Fermion)
    fatal_error("only fermionic Green's functions are supported");
   if(norm != 1.0) fatal_error("fermionic Green's functions must have a norm of 1.0");
   check_input_gf(g_tau,S_tau);
   if(!is_gf_real(g_tau) || !is_gf_real(S_tau))
    fatal_error("imaginary time Green's functions must be real");
   gf<imtime, matrix_real_valued> g_tau_real = real(g_tau), S_tau_real = real(S_tau);
   set_input_data(make_const_view(g_tau_real), make_const_view(S_tau_real));
  }
  break;
  case Susceptibility:
   // TODO
   fatal_error("continuation of susceptibilities is not yet implemented");
   break;
  case Conductivity:
   // TODO
   fatal_error("continuation of conductivity is not yet implemented");
   break;
  default:
   fatal_error("unknown observable kind " + to_string(kind));
 }
}

// Imaginary frequency
som_core::som_core(gf_const_view<imfreq> g_iw, gf_const_view<imfreq> S_iw,
                   observable_kind kind, double norm) :
 mesh(g_iw.mesh()), kind(kind), norm(norm),
 rhs(input_data_c_t()), error_bars(input_data_c_t()) {
 switch(kind) {
  case FermionGf: {
   if(g_iw.domain().statistic != Fermion)
    fatal_error("only fermionic Green's functions are supported");
   if(norm != 1.0) fatal_error("fermionic Green's functions must have a norm of 1.0");
   check_input_gf(g_iw,S_iw);
   if(!is_gf_real_in_tau(g_iw) || !is_gf_real_in_tau(S_iw))
    fatal_error("imaginary frequency Green's functions must correspond to a real G(\\tau)");
   auto g_iw_pos = positive_freq_view(g_iw);
   auto S_iw_pos = positive_freq_view(S_iw);
   check_input_gf(g_iw_pos,S_iw_pos);
   set_input_data(g_iw_pos,S_iw_pos);
  }
  break;
  case Susceptibility:
   // TODO
   fatal_error("continuation of susceptibilities is not yet implemented");
   break;
  case Conductivity:
   // TODO
   fatal_error("continuation of conductivity is not yet implemented");
   break;
  default:
   fatal_error("unknown observable kind " + to_string(kind));
 }
}

// Legendre coefficients
som_core::som_core(gf_const_view<legendre> g_l, gf_const_view<legendre> S_l,
                   observable_kind kind, double norm) :
 mesh(g_l.mesh()), kind(kind), norm(norm),
 rhs(input_data_r_t()), error_bars(input_data_r_t()) {
 switch(kind) {
  case FermionGf: {
   if(g_l.domain().statistic != Fermion)
    fatal_error("only fermionic Green's functions are supported");
   if(norm != 1.0) fatal_error("fermionic Green's functions must have a norm of 1.0");
   check_input_gf(g_l,S_l);
   if(!is_gf_real(g_l) || !is_gf_real(S_l))
    fatal_error("Legendre Green's functions must be real");
    gf<legendre, matrix_real_valued> g_l_real = real(g_l), S_l_real = real(S_l);
    set_input_data(make_const_view(g_l_real),make_const_view(S_l_real));
  }
  break;
  case Susceptibility:
   // TODO
   fatal_error("continuation of susceptibilities is not yet implemented");
   break;
  case Conductivity:
   // TODO
   fatal_error("continuation of conductivity is not yet implemented");
   break;
  default:
   fatal_error("unknown observable kind " + to_string(kind));
 }
}

void som_core::run(run_parameters_t const& p) {

 params = p;

 if((kind == Susceptibility || kind == Conductivity) && params.energy_window.first < 0) {
  params.energy_window.first = 0;
  if(params.verbosity > 0) std::cout << "WARNING: left boundary of the energy window is reset to 0";
 }

 if(params.energy_window.first >= params.energy_window.second)
  fatal_error("wrong energy window [" + to_string(params.energy_window.first) +
              ";" + to_string(params.energy_window.second) + "]");

 double window_width = params.energy_window.second - params.energy_window.first;
 if(params.min_rect_width <= 0 || params.min_rect_width >= window_width)
  fatal_error("min_rect_width must be in [0;" + to_string(window_width) + "]");

 if(params.min_rect_weight <= 0 || params.min_rect_weight >= norm)
  fatal_error("min_rect_weight must be in [0;" + to_string(norm) + "]");

 if(params.adjust_nsol_verygood_d > params.adjust_nsol_good_d)
  fatal_error("Cannot have adjust_nsol_verygood_d > adjust_nsol_good_d");

// FIXME: get rid of this macro
#define EI(ok, mk) int(ok) + 3 * mk

 switch(EI(kind, mesh.index())) {
  case EI(FermionGf,0): run_impl<kernel<FermionGf,imtime>>(); break;
  case EI(FermionGf,1): run_impl<kernel<FermionGf,imfreq>>(); break;
  //case EI(FermionGf,2): run_impl<kernel<FermionGf,legendre>>(); break;
  // TODO
 }
#undef EI
}

void som_core::operator()(gf_view<refreq> g_w) const {
 auto shape = get_target_shape(g_w);
 auto gf_dim = results.size();

 if(shape[0] != gf_dim || shape[1] != gf_dim)
  fatal_error("expected a real-frequency Green's function with matrix dimensions "
              + to_string(gf_dim) + "x" + to_string(gf_dim) + " in assignment");

 g_w() = 0;
 for(int i = 0; i < gf_dim; ++i) {
  auto const& conf = results[i];
  for(auto const& rect : conf) {
   for(auto e : g_w.mesh()) g_w.data()(e.index(),i,i) += rect.hilbert_transform(double(e));
   auto & tail = g_w.singularity();
   tail.data()(range(),i,i) += rect.tail_coefficients(tail.order_min(),tail.order_max());
  }
 }
}

template<typename KernelType> void som_core::run_impl() {

 using mesh_t = typename KernelType::mesh_type;
 mesh_t const& m = mesh;

 if(params.verbosity > 0) {
  std::cout << "Constructing integral kernel... ";
 }
 KernelType kernel(m);
 if(params.verbosity > 0) {
  std::cout << "done" << std::endl;
  std::cout << "Kernel: " << kernel << std::endl;
 }

 // Find solution for each component of GF
 for(int i = 0; i < results.size(); ++i) {
  if(params.verbosity > 0)
   std::cout << "Running algorithm for observable component [" << i << "," << i << "]" << std::endl;

  auto const& rhs_ = (input_data_t<mesh_t> const&) rhs;
  auto const& error_bars_ = (input_data_t<mesh_t> const&) error_bars;

  int F = params.adjust_ngu ? adjust_f(kernel,rhs_[i],error_bars_[i]) : params.n_global_updates;
  results[i] = accumulate(kernel,rhs_[i],error_bars_[i],histograms[i],F);
 }
}

template<typename KernelType> int som_core::adjust_f(KernelType const& kern,
 typename KernelType::result_type rhs_,
 typename KernelType::result_type error_bars_) {

 if(params.verbosity >= 1) {
  std::cout << "Adjusting the number of global updates F using "
            << params.adjust_ngu_n_solutions << " particular solutions ..." << std::endl;
 }
 objective_function<KernelType> of(kern, rhs_, error_bars_);
 fit_quality<KernelType> fq(kern, rhs_, error_bars_);

 int F = params.n_global_updates;

 int n_good_solutions;
 for(n_good_solutions = 0;;n_good_solutions = 0, F *= 2) {
  solution_worker<KernelType> worker(of,norm,ci,params,F);
  auto & rng = worker.get_rng();

  int n_sol;
  for(int i = 0; (n_sol = comm.rank() + i*comm.size()) < params.adjust_ngu_n_solutions; ++i) {
   if(params.verbosity >= 2) {
    std::cout << "[Node " << comm.rank() << "] Accumulation of particular solution "
              << n_sol << std::endl;
   }

   auto solution = worker(1 + rng(params.max_rects));
   double kappa = fq(solution);

   if(kappa > params.adjust_ngu_kappa) ++n_good_solutions;
   if(params.verbosity >= 2) {
    std::cout << "[Node " << comm.rank() << "] Particular solution "
              << n_sol << " is " << (kappa > params.adjust_ngu_kappa ? "" : "not ")
              << "good (\\kappa = " << kappa
              << ", D = " << worker.get_objf_value() << ")." << std::endl;
   }
  }
  comm.barrier();
  n_good_solutions = mpi_all_reduce(n_good_solutions,comm);

  if(params.verbosity >= 1)
   std::cout << "F = " << F << ", "
             << n_good_solutions << " solutions with \\kappa > " << params.adjust_ngu_kappa
             << " (out of " << params.adjust_ngu_n_solutions << ")" << std::endl;

  if(n_good_solutions > params.adjust_ngu_n_solutions/2) break;
 }
 if(params.verbosity >= 1) std::cout << "F = " << F << " is enough." << std::endl;

 return F;
}

template<typename KernelType> configuration som_core::accumulate(KernelType const& kern,
 typename KernelType::result_type rhs_,
 typename KernelType::result_type error_bars_,
 histogram & hist, int F) {

 if(params.verbosity >= 1) std::cout << "Accumulating particular solutions ..." << std::endl;

 objective_function<KernelType> of(kern, rhs_, error_bars_);
 solution_worker<KernelType> worker(of,norm,ci,params,F);
 auto & rng = worker.get_rng();

 // Pairs (configuration, objective function)
 std::vector<std::pair<configuration,double>> solutions;

 int n_sol_max = 0;                          // Number of solutions to be accumulated
 int n_sol, i = 0;                           // Global and rank-local indices of solution
 int n_good_solutions, n_verygood_solutions; // Number of good and very good solutions
 double objf_min = HUGE_VAL;                 // Minimal value of D
 do {
  n_sol_max += params.n_solutions;
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
   if(s.second/objf_min <= params.adjust_nsol_good_d) ++n_good_solutions;
   if(s.second/objf_min <= params.adjust_nsol_verygood_d) ++n_verygood_solutions;
  }
  n_good_solutions = mpi_all_reduce(n_good_solutions);
  n_verygood_solutions = mpi_all_reduce(n_verygood_solutions);

  if(params.verbosity >= 1) {
   std::cout << "D_min = " << objf_min << std::endl;
   std::cout << "Number of good soulutions (D/D_min <= "
             << params.adjust_nsol_good_d << ") = " << n_good_solutions << std::endl;
   std::cout << "Number of very good soulutions (D/D_min <= "
             << params.adjust_nsol_verygood_d << ") = " << n_verygood_solutions << std::endl;
  }

 } while(params.adjust_nsol &&
         double(n_verygood_solutions) / double(n_good_solutions) < params.adjust_nsol_ratio);

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
  if(s.second/objf_min <= params.adjust_nsol_good_d) sol_sum += s.first;
 }
 sol_sum *= double(1/n_good_solutions);

 // Sum over all processes
 sol_sum = mpi_reduce(sol_sum, comm, 0, true);

 if(params.make_histograms) hist = mpi_reduce(hist, comm, 0, true);

 if(params.verbosity >= 1) std::cout << "Done" << std::endl;

 return sol_sum;
}

som_core::~som_core() {
 int n_configs = results.size();
 int n_rects = 0;
 for(auto const& c : results) n_rects += c.size();
 if(n_configs + n_rects != ci.n_used_entries()) TRIQS_RUNTIME_ERROR <<
  "Detected inconsistency in LHS cache: " <<
  "n_configs = " << n_configs << ", n_rects = " << n_rects <<
  ", n_used_entries = " << ci.n_used_entries();
}

}
