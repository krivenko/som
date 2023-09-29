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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>

#include <itertools/itertools.hpp>

#include <nda/linalg.hpp>

#include <som/global_index_map.hpp>
#include <som/kernels/all.hpp>
#include <som/numerics/ecqp_worker.hpp>
#include <som/numerics/finite_diff.hpp>
#include <som/solution_functionals/objective_function.hpp>

#include "common.hxx"
#include "som_core.hpp"

namespace som {

using namespace nda;
using namespace itertools;

////////////////////////////////////////////
// som_core::data_t::compute_objf_final() //
////////////////////////////////////////////

template <typename KernelType>
void som_core::data_t::compute_objf_final(KernelType const& kernel) {
  auto of = make_objf<KernelType>(kernel);
  objf_final = of(final_solution);
}

////////////////////////////////////
// som_core::compute_objf_final() //
////////////////////////////////////

template <typename KernelType>
std::vector<double> som_core::compute_objf_final_impl() {
  using mesh_t = typename KernelType::mesh_type;
  mesh_t const& m = std::get<mesh_t>(mesh);

  KernelType kernel(m);

  std::vector<double> objf_final(data.size());
  for(auto n : range(long(data.size()))) {
    auto& d = data[n];
    d.compute_objf_final(kernel);
    objf_final[n] = d.objf_final;
  }

  ci.invalidate_all();

  return objf_final;
}

std::vector<double> som_core::compute_objf_final() {
#define IMPL_CASE(r, okmk)                                                     \
  case(kernel_id<BOOST_PP_SEQ_ELEM(1, okmk)>(BOOST_PP_SEQ_ELEM(0, okmk))):     \
    return compute_objf_final_impl<kernel<BOOST_PP_SEQ_ENUM(okmk)>>();

  SELECT_KERNEL(IMPL_CASE, som_core::compute_objf_final())
#undef IMPL_CASE
}

////////////////////////////////////////
// som_core::compute_final_solution() //
////////////////////////////////////////

std::vector<double> som_core::compute_final_solution(double good_chi_rel,
                                                     double good_chi_abs,
                                                     int verbosity) {
  for(auto n : range(long(data.size()))) {
    if(verbosity > 0) {
      mpi_cout(comm) << "Constructing the final solution "
                        "for observable component ["
                     << n << "," << n << "]" << std::endl;
    }

    auto& d = data[n];

    configuration sol_sum(ci);

    // Rank-local stage of summation
    std::size_t n_good_solutions = 0;
    double const chi_min = std::sqrt(d.objf_min);
    double const chi_c = std::min(good_chi_abs, chi_min * good_chi_rel);
    if(verbosity > 0) {
      mpi_cout(comm) << "Good solution threshold χ_c = " << chi_c << std::endl;
    }

    for(auto const& s : d.particular_solutions) {
      // Pick only good solutions
      double chi = std::sqrt(s.second);
      if(chi <= chi_c) {
        sol_sum += s.first;
        ++n_good_solutions;
      }
    }

    // Sum over all processes
    n_good_solutions = mpi::all_reduce(n_good_solutions, comm);

    if(n_good_solutions == 0)
      TRIQS_RUNTIME_ERROR
          << "No good solution could be selected, try accumulating more "
             "solutions, and/or set a higher threshold χ_c by increasing "
             "values of good_chi_rel/good_chi_abs";
    else {
      if(verbosity > 0) {
        mpi_cout(comm) << "Selected " << n_good_solutions << " good solutions"
                       << std::endl;
      }
    }

    d.final_solution = mpi::all_reduce(sol_sum, comm);
    d.final_solution *= 1.0 / double(n_good_solutions);
  }

  return compute_objf_final();
}

///////////////////////////////////////////
// som_core::compute_final_solution_cc() //
///////////////////////////////////////////

void final_solution_cc_parameters_t::validate() const {
  using std::to_string;

  struct {
    std::size_t mesh_size = {};
    void operator()(triqs::mesh::refreq const& m) {
      mesh_size = m.size();
      if(mesh_size < 3)
        fatal_error("There must be at least 3 energy points in refreq_mesh");
    }
    void operator()(array<double, 1> const& m) {
      mesh_size = m.size();
      if(mesh_size < 3)
        fatal_error("There must be at least 3 energy points in refreq_mesh");
      if(!std::is_sorted(m.begin(), m.end(), std::less_equal<>()))
        fatal_error(
            "Energy points in refreq_mesh must be listed in strictly ascending "
            "order");
    }
  } validate_mesh;
  std::visit(validate_mesh, refreq_mesh);

  if(good_chi_rel < 1.0)
    fatal_error("Parameter good_chi_rel must be >= 1.0 (got " +
                to_string(good_chi_rel) + ")");

  if(good_chi_abs < 0)
    fatal_error("Parameter good_chi_abs must be non-negative (got " +
                to_string(good_chi_abs) + ")");

  if((!default_model.empty()) &&
     (default_model.size() != validate_mesh.mesh_size))
    fatal_error("default_model must have the same size as refreq_mesh (got " +
                to_string(default_model.size()) + " vs " +
                to_string(validate_mesh.mesh_size) + ")");

  if((!default_model_weights.empty()) &&
     (default_model_weights.size() != validate_mesh.mesh_size))
    fatal_error(
        "default_model_weights must have the same size as refreq_mesh (got " +
        to_string(default_model_weights.size()) + " vs " +
        to_string(validate_mesh.mesh_size) + ")");

  if(max_iter <= 0)
    fatal_error("Number of iterations max_iter must be positive (got " +
                to_string(max_iter) + ")");

  if(ew_penalty_coeff < 0)
    fatal_error("Parameter ew_penalty_coeff must positive (got " +
                to_string(ew_penalty_coeff) + ")");

  if(amp_penalty_max <= 0)
    fatal_error("Parameter amp_penalty_max must be positive (got " +
                to_string(amp_penalty_max) + ")");

  if(amp_penalty_divisor <= 1.0)
    fatal_error("Parameter amp_penalty_divisor must exceed 1.0 (got " +
                to_string(amp_penalty_divisor) + ")");

  if(der_penalty_init <= 0)
    fatal_error("Parameter der_penalty_init must be positive (got " +
                to_string(der_penalty_init) + ")");

  if(der_penalty_coeff <= 1.0)
    fatal_error("Parameter der_penalty_coeff must exceed 1.0 (got " +
                to_string(der_penalty_coeff) + ")");
}

vector<double> make_rhs_vector(array<double, 2> const& A_j_block,
                               double U,
                               array<double, 1> const& T_k,
                               array<double, 1> const& A_T_k,
                               global_index_map const& index_map,
                               mpi::communicator const& comm) {
  auto const J = index_map.size();
  auto f = vector<double>::zeros({long(J)});

  // RHS generated by O_T
  if(!T_k.empty()) {
    int const my_rank = comm.rank();
    auto const my_j_start = long(index_map.range_start(my_rank));
    auto const my_j_size = long(index_map.range_size(my_rank));
    for(auto j : range(my_j_size))
      f(my_j_start + j) += sum(T_k * A_T_k * A_j_block(j, range::all));
    f = mpi::all_reduce(f, comm);
  }

  // RHS generated by O_U
  f() += U / double(J);

  return f;
}

void update_O_mat(array<double, 2>& A_j_local_block,
                  array<double, 2>& Ap_j_local_block,
                  array<double, 2>& App_j_local_block,
                  double U,
                  array<double, 1> const& Q_k,
                  array<double, 1> const& D_k,
                  array<double, 1> const& B_k,
                  array<double, 1> const& T_k,
                  global_index_map const& index_map,
                  mpi::communicator const& comm,
                  bool verbose,
                  matrix<double, nda::F_layout>& O_mat) {
  int const j_block1 = comm.rank();
  auto const J1 = A_j_local_block.shape()[0];
  assert(Ap_j_local_block.shape()[0] == J1);
  assert(App_j_local_block.shape()[0] == J1);

  if(verbose)
    mpi_cout(comm) << "    Updating matrix of quadratic form O" << std::endl;

  // A_j(e_k) block received from other ranks
  array<double, 2> A_j_remote_block;
  // A'_j(e_k) block received from other ranks
  array<double, 2> Ap_j_remote_block;
  // A''_j(e_k) block received from other ranks
  array<double, 2> App_j_remote_block;

  auto const& A_j_block1 = A_j_local_block;
  auto const& Ap_j_block1 = Ap_j_local_block;
  auto const& App_j_block1 = App_j_local_block;

  O_mat() = 0;

  // NOLINTNEXTLINE(cppcoreguidelines-narrowing-conversions)
  for(int j_block2 : range(comm.size())) {
    bool diag_block = j_block2 == j_block1;

    auto& A_j_block2 = diag_block ? A_j_local_block : A_j_remote_block;
    auto& Ap_j_block2 = diag_block ? Ap_j_local_block : Ap_j_remote_block;
    auto& App_j_block2 = diag_block ? App_j_local_block : App_j_remote_block;

    mpi::broadcast(A_j_block2, comm, j_block2);
    mpi::broadcast(Ap_j_block2, comm, j_block2);
    mpi::broadcast(App_j_block2, comm, j_block2);

    auto const J2 = A_j_block2.shape()[0];
    assert(Ap_j_block2.shape()[0] == J2);
    assert(App_j_block2.shape()[0] == J2);

    if(verbose) {
      mpi_cout(comm) << "    Combining " << J1 << " solutions from rank "
                     << j_block1 << " with " << J2 << " solutions from rank "
                     << j_block2 << std::endl;
    }

    for(auto j1_local : range(J1)) {
      auto j1 = index_map(j_block1, j1_local);
      for(auto j2_local : range(J2)) {
        auto j2 = index_map(j_block2, j2_local);

        // Terms O_Q, O_D and O_B
        O_mat(j1, j2) += sum(Q_k * Q_k * A_j_block1(j1_local, range::all) *
                             A_j_block2(j2_local, range::all));
        O_mat(j1, j2) += sum(D_k * D_k * Ap_j_block1(j1_local, range::all) *
                             Ap_j_block2(j2_local, range::all));
        O_mat(j1, j2) += sum(B_k * B_k * App_j_block1(j1_local, range::all) *
                             App_j_block2(j2_local, range::all));

        // Term O_T
        if(!T_k.empty())
          O_mat(j1, j2) += sum(T_k * A_j_block1(j1_local, range::all) *
                               A_j_block2(j2_local, range::all));
      }
    }
  }

  O_mat = mpi::all_reduce(O_mat, comm);

  // Penalty for large deviations from the equal-weight superposition, O_U
  O_mat += U;
}

struct cc_protocol_iter_result_t {
  int iter;                 // Iteration number
  nda::vector<double> c;    // Coefficients c_j
  double sum_abs;           // \sum_j |c_j|
  double diff_c_normalized; // \sum_j |c^iter_j - c^{iter-1}_j| / \sum_j |c_j|
};

std::vector<cc_protocol_iter_result_t> cc_protocol(
    std::vector<std::pair<configuration, double>> const& particular_solutions,
    std::vector<std::size_t> const& good_solution_indices,
    array<double, 1> const e_k_list,
    double convergence_tol,
    final_solution_cc_parameters_t const& p,
    global_index_map const& index_map,
    mpi::communicator const& comm) {
  auto const my_J = good_solution_indices.size();

  int const my_rank = comm.rank();
  auto const my_j_start = long(index_map.range_start(my_rank));
  auto const my_j_size = long(index_map.range_size(my_rank));
  auto const my_j_range = range(my_j_size);

  // Rank-local portion of A_j(e_k) array (partition over the first index)
  array<double, 2> A_block(my_J, e_k_list.size());
  // Rank-local portion of A'_j(e_k) array (partition over the first index)
  array<double, 2> Ap_block(my_J, e_k_list.size() - 1);
  // Rank-local portion of A''_j(e_k) array (partition over the first index)
  array<double, 2> App_block(my_J, e_k_list.size() - 2);

  for(auto const& [j, s_index] : enumerate(good_solution_indices)) {
    // Fill A
    for(auto [k, e_k] : enumerate(e_k_list))
      A_block(j, k) = particular_solutions[s_index].first(e_k);

    auto A_view = make_array_const_view(A_block(j, range::all));
    // Fill Ap
    finite_diff_forward(A_view, e_k_list, Ap_block(j, range::all));
    // Fill App
    finite_diff_2_symm(A_view, e_k_list, App_block(j, range::all));
  }

  // Sum of A_j(e_k) over j
  array<double, 1> A_k(A_block.shape()[1]);
  // Sum of A'_j(e_k) over j
  array<double, 1> Ap_k(Ap_block.shape()[1]);
  // Sum of A''_j(e_k) over j
  array<double, 1> App_k(App_block.shape()[1]);

  // Total number of selected particular solutions (summed over all ranks)
  auto const J = index_map.size();

  double DD = p.der_penalty_init / double(J);

  // Amplitude regularization parameters Q_k
  array<double, 1> Q_k = zeros<double>(A_k.size());
  // First derivative regularization parameters D_k
  array<double, 1> D_k = DD * ones<double>(Ap_k.size());
  // Second derivative regularization parameters B_k
  array<double, 1> B_k = DD * ones<double>(App_k.size());

  // List of coefficients c_j, one element per iteration
  std::vector<cc_protocol_iter_result_t> results = {
      {0, (1.0 / double(J)) * ones<double>(J), 1.0, HUGE_VAL}};

  // Prepare RHS vector
  auto f = make_rhs_vector(A_block,
                           p.ew_penalty_coeff,
                           p.default_model_weights,
                           p.default_model,
                           index_map,
                           comm);

  // Matrix of the quadratic form
  matrix<double, nda::F_layout> O_mat(J, J);

  // Equality Constrained Quadratic Programming worker
  nda::matrix<double, nda::F_layout> c_sum_cons_matrix =
      nda::ones<double>(1, J);
  nda::vector<double> c_sum_cons_rhs = {1};
  auto worker = ecqp_worker(int(J), 1);

  int iter = 0;
  for(; iter < p.max_iter; ++iter) {
    if(p.verbosity >= 1)
      mpi_cout(comm) << "  Iteration " << (iter + 1) << " out of " << p.max_iter
                     << std::endl;

    // Update matrix of the quadratic form
    update_O_mat(A_block,
                 Ap_block,
                 App_block,
                 p.ew_penalty_coeff,
                 Q_k,
                 D_k,
                 B_k,
                 p.default_model_weights,
                 index_map,
                 comm,
                 p.verbosity >= 2,
                 O_mat);

    results.emplace_back(
        cc_protocol_iter_result_t{iter + 1, nda::vector<double>(J), 0.0, 0.0});

    auto& res = results.back();
    auto const my_c_view = res.c(range(my_j_start, my_j_start + my_j_size));

    // Minimize the functional
    worker(O_mat(), f, c_sum_cons_matrix(), c_sum_cons_rhs(), res.c);

    // Linear combination of particular solutions: Rank-local stage
    if(my_J != 0) {
      nda::clef::placeholder<0> j_;
      nda::clef::placeholder<1> k_;
      using nda::clef::sum;
      A_k(k_) << sum(my_c_view(j_) * A_block(j_, k_), j_ = my_j_range);
      Ap_k(k_) << sum(my_c_view(j_) * Ap_block(j_, k_), j_ = my_j_range);
      App_k(k_) << sum(my_c_view(j_) * App_block(j_, k_), j_ = my_j_range);
    } else {
      A_k() = 0;
      Ap_k() = 0;
      App_k() = 0;
    }

    // Linear combination of particular solutions: Inter-rank stage
    A_k = mpi::all_reduce(A_k, comm);
    Ap_k = mpi::all_reduce(Ap_k, comm);
    App_k = mpi::all_reduce(App_k, comm);

    res.sum_abs = sum(abs(res.c));
    res.diff_c_normalized = sum(abs(res.c - results[iter].c)) / res.sum_abs;

    std::string sum_c_text, diff_c_text;
    if(p.verbosity >= 1) {
      sum_c_text = "sum(|c_" + std::to_string(iter + 1) + "|)";
      diff_c_text = "sum(|c_" + std::to_string(iter + 1) + " - c_" +
                    std::to_string(iter) + "|)";
    }

    if(p.monitor) {
      bool stop = p.monitor(res.c, {A_k, Q_k}, {Ap_k, D_k}, {App_k, B_k});
      if(stop) {
        if(p.verbosity >= 1) {
          mpi_cout(comm)
              << "  Monitor function returned true, stopping iterations: "
              << sum_c_text << " = " << res.sum_abs << ", " << diff_c_text
              << " / " << sum_c_text << " = " << res.diff_c_normalized
              << std::endl;
        }
        break;
      }
    }

    if(p.verbosity >= 1) {
      mpi_cout(comm) << "  End of iteration " << (iter + 1) << ": "
                     << sum_c_text << " = " << res.sum_abs << ", "
                     << diff_c_text << " / " << sum_c_text << " = "
                     << res.diff_c_normalized << std::endl;
    }

    if(res.sum_abs > p.max_sum_abs_c) {
      if(p.verbosity >= 1)
        mpi_cout(comm)
            << "Failed to converge, " << sum_c_text << " = " << res.sum_abs
            << " is too big. "
            << "Consider decreasing der_penalty_coeff and/or der_penalty_init."
            << std::endl;
      break;
    }

    if(iter > 0 && (res.diff_c_normalized < convergence_tol)) {
      if(p.verbosity >= 1) mpi_cout(comm) << "Convergence reached" << std::endl;
      break;
    }

    //
    // Update regularization parameters
    //

    DD *= p.der_penalty_coeff;
    D_k *= p.der_penalty_coeff;
    B_k *= p.der_penalty_coeff;

    // Update Q
    for(auto k : range(A_k.size())) {
      if(A_k(k) < 0)
        Q_k(k) = p.amp_penalty_max;
      else
        Q_k(k) /= p.amp_penalty_divisor;
    }

    // Update D
    for(auto k : range(Ap_k.size())) {
      double val = std::abs(Ap_k(k));
      if(D_k(k) * val > DD) D_k(k) = DD / val;
    }

    // Update B
    for(auto k : range(App_k.size())) {
      double val = std::abs(App_k(k));
      if(B_k(k) * val > DD) B_k(k) = DD / val;
    }
  }

  if(p.verbosity >= 1 && iter == p.max_iter)
    mpi_cout(comm) << "Maximum number of iterations reached" << std::endl;

  return results;
}

template <typename KernelType>
std::vector<double> som_core::compute_final_solution_cc_impl(
    final_solution_cc_parameters_t const& p) {

  // List of real frequency points
  auto e_k_list = std::visit(
      [](auto const& m) {
        array<double, 1> res(m.size());
        for(auto const& [k, e] : enumerate(m)) res(k) = double(e);
        return res;
      },
      p.refreq_mesh);

  for(auto n : range(long(data.size()))) {
    auto& d = data[n];

    if(p.verbosity > 0)
      mpi_cout(comm) << "Using CC protocol to construct the final solution "
                        "for observable component ["
                     << n << "," << n << "]" << std::endl;

    // Select good solutions
    std::vector<std::size_t> good_solution_indices;
    good_solution_indices.reserve(d.particular_solutions.size());
    double const chi_min = std::sqrt(d.objf_min);
    double const chi_c = std::min(p.good_chi_abs, chi_min * p.good_chi_rel);

    if(p.verbosity > 0) {
      mpi_cout(comm) << "Good solution threshold χ_c = " << chi_c << std::endl;
    }

    auto is_good = [chi_c](double chi2) { return std::sqrt(chi2) <= chi_c; };

    for(auto const& [s_index, s] : enumerate(d.particular_solutions)) {
      if(is_good(s.second)) good_solution_indices.emplace_back(s_index);
    }

    global_index_map index_map(comm, good_solution_indices.size());
    if(index_map.size() == 0) {
      TRIQS_RUNTIME_ERROR
          << "No good solution could be selected, try accumulating more "
             "solutions, and/or set a higher threshold χ_c by increasing "
             "values "
             "of good_chi_rel/good_chi_abs";
    } else {
      if(p.verbosity > 0)
        mpi_cout(comm) << "Selected " << good_solution_indices.size()
                       << " good solutions on this rank and "
                       << index_map.size() << " good solutions in total"
                       << std::endl;
    }

    using mesh_t = typename KernelType::mesh_type;
    auto const& m = std::get<mesh_t>(mesh);
    KernelType kernel(m);
    auto of = d.make_objf<KernelType>(kernel);

    // Convergence tolerance for CC iterations
    auto const& U_dagger = of.get_U_dagger();
    auto rhs =
        U_dagger
            ? (*U_dagger) *
                  vector_const_view<typename decltype(of)::rhs_scalar_type>(
                      of.get_rhs())
            : of.get_rhs();
    double convergence_tol = min_element(sqrt(of.get_sigma2()) / abs(rhs));

    if(p.verbosity > 1)
      mpi_cout(comm) << "Convergence tolerance is " << convergence_tol
                     << std::endl;

    // Compute coefficients c_j using CC optimization protocol
    auto results = cc_protocol(d.particular_solutions,
                               good_solution_indices,
                               e_k_list,
                               convergence_tol,
                               p,
                               index_map,
                               comm);

    configuration sol(ci);
    auto const j_range_start = index_map.range_start(comm.rank());

    // Order elements of 'results' according to diff_c_normalized, e.i.
    // to proximity to convergence.
    std::sort(
        results.begin(), results.end(), [](auto const& r1, auto const& r2) {
          return r1.diff_c_normalized < r2.diff_c_normalized;
        });

    // Pick the first good final solution from sorted 'results'
    int selected_iter = -1;
    const nda::vector<double>* selected_c = nullptr;
    for(auto const& res : results) {
      for(auto const& [j, s_index] : enumerate(good_solution_indices)) {
        sol += res.c(j_range_start + j) * d.particular_solutions[s_index].first;
      }
      sol = mpi::all_reduce(sol, comm);

      double chi2 = of(sol);
      if(p.verbosity >= 2) {
        mpi_cout(comm)
            << "Linear combination of particular solutions from iteration "
            << res.iter << ": χ = " << std::sqrt(chi2) << std::endl;
      }
      if(is_good(chi2)) {
        selected_iter = res.iter;
        selected_c = &res.c;
        break;
      }

      sol.clear();
    }

    if(selected_iter == -1) {
      ci.invalidate_all();
      TRIQS_RUNTIME_ERROR << "Could not construct the final solution, try "
                             "increasing values of good_chi_rel/good_chi_abs";
    }

    if(p.verbosity >= 1) {
      mpi_cout(comm)
          << "Forming the resulting linear combination of particular solutions "
             "using coefficients from iteration "
          << selected_iter << std::endl;
      if(p.verbosity >= 2)
        mpi_cout(comm)
            << "Coefficients of the linear combination from iteration "
            << selected_iter << ": " << *selected_c << std::endl;
    }

    d.final_solution = sol;
  }

  ci.invalidate_all();

  return compute_objf_final();
}

std::vector<double>
som_core::compute_final_solution_cc(final_solution_cc_parameters_t const& p) {
  p.validate();

#define IMPL_CASE(r, okmk)                                                     \
  case(kernel_id<BOOST_PP_SEQ_ELEM(1, okmk)>(BOOST_PP_SEQ_ELEM(0, okmk))):     \
    return compute_final_solution_cc_impl<kernel<BOOST_PP_SEQ_ENUM(okmk)>>(p);

  SELECT_KERNEL(IMPL_CASE, som_core::compute_final_solution_cc())
#undef IMPL_CASE
}

/////////////////////
// som_core::run() //
/////////////////////

void som_core::run(accumulate_parameters_t const& p) {
  accumulate(p);
  compute_final_solution(p.adjust_l_good_chi);
}

} // namespace som
