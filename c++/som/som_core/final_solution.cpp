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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

#include <itertools/itertools.hpp>

#include <triqs/arrays/linalg/det_and_inverse.hpp>
#include <triqs/clef.hpp>

#include <som/global_index_map.hpp>
#include <som/kernels/all.hpp>
#include <som/numerics/finite_diff.hpp>
#include <som/solution_functionals/objective_function.hpp>

#include "common.hxx"
#include "som_core.hpp"

namespace som {

using namespace triqs::arrays;
using namespace itertools;

////////////////////////////////////////////
// som_core::data_t::compute_objf_final() //
////////////////////////////////////////////

template <typename KernelType>
void som_core::data_t::compute_objf_final(KernelType const& kernel) {
  using mesh_t = typename KernelType::mesh_type;
  objective_function<KernelType> of(
      kernel, get_rhs<mesh_t>(), get_error_bars<mesh_t>());
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
  for(int n : range(data.size())) {
    auto& d = data[n];
    d.compute_objf_final(kernel);
    objf_final[n] = d.objf_final;
  }

  ci.invalidate_all();

  return objf_final;
}

std::vector<double> som_core::compute_objf_final() {
#define IMPL_CASE(r, okmk)                                                     \
  case(kernel_id(BOOST_PP_SEQ_ELEM(0, okmk), BOOST_PP_SEQ_ELEM(1, okmk){})):   \
    return compute_objf_final_impl<kernel<BOOST_PP_SEQ_ENUM(okmk)>>();

  SELECT_KERNEL(IMPL_CASE, som_core::compute_objf_final())
#undef IMPL_CASE
}

////////////////////////////////////////
// som_core::compute_final_solution() //
////////////////////////////////////////

std::vector<double> som_core::compute_final_solution(double good_chi_rel,
                                                     double good_chi_abs) {
  for(int n : range(data.size())) {
    auto& d = data[n];

    configuration sol_sum(ci);

    // Rank-local stage of summation
    int n_good_solutions = 0;
    double chi_min = std::sqrt(d.objf_min);
    for(auto const& s : d.particular_solutions) {
      // Pick only good solutions
      double chi = std::sqrt(s.second);
      if(chi / chi_min <= good_chi_rel && chi <= good_chi_abs) {
        sol_sum += s.first;
        ++n_good_solutions;
      }
    }

    // Sum over all processes
    n_good_solutions = mpi::all_reduce(n_good_solutions, comm);
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
    int mesh_size;
    void operator()(triqs::gfs::gf_mesh<triqs::gfs::refreq> const& m) {
      mesh_size = m.size();
      if(mesh_size < 3)
        fatal_error("There must be at least 3 energy points in refreq_mesh");
    }
    void operator()(array<double, 1> const& m) {
      mesh_size = m.size();
      if(mesh_size < 3)
        fatal_error("There must be at least 3 energy points in refreq_mesh");
      if(!std::is_sorted(m.begin(), m.end(), std::less_equal<double>()))
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

  if((default_model.size() != 0) && // TODO: Use nda::array::empty()
     (default_model.size() != validate_mesh.mesh_size))
    fatal_error("default_model must have the same size as refreq_mesh (got " +
                to_string(default_model.size()) + " vs " +
                to_string(validate_mesh.mesh_size) + ")");

  if((default_model_weights.size() != 0) && // TODO: Use nda::array::empty()
     (default_model_weights.size() != validate_mesh.mesh_size))
    fatal_error(
        "default_model_weights must have the same size as refreq_mesh (got " +
        to_string(default_model_weights.size()) + " vs " +
        to_string(validate_mesh.mesh_size) + ")");

  if(max_iter <= 0)
    fatal_error("Number of iterations max_iter must be positive (got " +
                to_string(max_iter) + ")");

  if(unity_sum_coeff < 0)
    fatal_error("Parameter unity_sum_coeff must positive (got " +
                to_string(unity_sum_coeff) + ")");

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
  int J = index_map.size();
  vector<double> b(J);
  b() = 0;

  // RHS generated by O_5
  if(T_k.size() != 0) { // TODO: Use nda::array::empty()
    for(int j : range(first_dim(A_j_block)))
      b(j) += sum(T_k * A_T_k * A_j_block(j, range()));
    b = mpi::all_reduce(b, comm);
  }

  // RHS generated by O_2 and O_3
  b() += U + 1.0 / double(J);

  return b;
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
                  matrix<double>& O_mat) {
  int const j_block1 = comm.rank();
  int const J1 = first_dim(A_j_local_block);
  assert(first_dim(Ap_j_local_block) == J1);
  assert(first_dim(App_j_local_block) == J1);

  if(verbose)
    mpi_cout(comm) << "    Updating matrix of quadratic form O" << std::endl;

  // Terms O_1, O_4 and O_5

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

  for(int j_block2 : range(comm.size())) {
    bool diag_block = j_block2 == j_block1;

    auto& A_j_block2 = diag_block ? A_j_local_block : A_j_remote_block;
    auto& Ap_j_block2 = diag_block ? Ap_j_local_block : Ap_j_remote_block;
    auto& App_j_block2 = diag_block ? App_j_local_block : App_j_remote_block;

    mpi::broadcast(A_j_block2, comm, j_block2);
    mpi::broadcast(Ap_j_block2, comm, j_block2);
    mpi::broadcast(App_j_block2, comm, j_block2);

    int const J2 = first_dim(A_j_block2);
    assert(first_dim(Ap_j_block2) == J2);
    assert(first_dim(App_j_block2) == J2);

    if(verbose) {
      mpi_cout(comm) << "    Combining " << J1 << " solutions from rank "
                     << j_block1 << " with " << J2 << " solutions from rank "
                     << j_block2 << std::endl;
    }

    for(int j1_local : range(J1)) {
      int j1 = index_map(j_block1, j1_local);
      for(int j2_local : range(J2)) {
        int j2 = index_map(j_block2, j2_local);

        O_mat(j1, j2) += sum(Q_k * A_j_block1(j1_local, range()) *
                             A_j_block2(j2_local, range()));
        O_mat(j1, j2) += sum(D_k * D_k * Ap_j_block1(j1_local, range()) *
                             Ap_j_block2(j2_local, range()));
        O_mat(j1, j2) += sum(B_k * B_k * App_j_block1(j1_local, range()) *
                             App_j_block2(j2_local, range()));

        if(T_k.size() != 0) // TODO: Use nda::array::empty()
          O_mat(j1, j2) += sum(T_k * A_j_block1(j1_local, range()) *
                               A_j_block2(j2_local, range()));
      }
    }
  }

  O_mat = mpi::all_reduce(O_mat, comm);

  // Unity sum constraint O_2
  O_mat.as_array_view() += U;

  // Penalty for large deviations from the equal-weight superposition, O_3
  O_mat() += 1.0;
}

vector<double> cc_protocol(
    std::vector<std::pair<configuration, double>> const& particular_solutions,
    std::vector<int> const& good_solution_indices,
    array<double, 1> const e_k_list,
    double convergence_tol,
    final_solution_cc_parameters_t const& p,
    global_index_map const& index_map,
    mpi::communicator const& comm) {
  int my_J = good_solution_indices.size();

  int const my_rank = comm.rank();
  auto my_j_start = index_map.range_start(my_rank);
  auto my_j_size = index_map.range_size(my_rank);

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
    auto A_view = make_const_view(A_block(j, range()));

    // Fill Ap
    finite_diff_forward(A_view, e_k_list, Ap_block(j, range()));
    // Fill App
    finite_diff_2_symm(A_view, e_k_list, App_block(j, range()));
  }

  // Sum of A_j(e_k) over j
  array<double, 1> A_k(second_dim(A_block));
  // Sum of A'_j(e_k) over j
  array<double, 1> Ap_k(second_dim(Ap_block));
  // Sum of A''_j(e_k) over j
  array<double, 1> App_k(second_dim(App_block));

  double DD = p.der_penalty_init;

  // Amplitude regularization parameters Q_k
  array<double, 1> Q_k(A_k.size());
  Q_k() = 0; // TODO: Use nda::zeros()
  // First derivative regularization parameters D_k
  array<double, 1> D_k(Ap_k.size());
  D_k() = DD; // TODO: Use nda::ones()
  // Second derivative regularization parameters B_k
  array<double, 1> B_k(App_k.size());
  B_k() = DD; // TODO: Use nda::ones()

  // Total number of selected particular solutions (summed over all ranks)
  int J = index_map.size();

  // Coefficients c_j
  vector<double> c(J);
  c() = 1.0 / double(J); // TODO: Use nda::ones()

  // Prepare RHS vector
  auto b = make_rhs_vector(A_block,
                           p.unity_sum_coeff,
                           p.default_model_weights,
                           p.default_model,
                           index_map,
                           comm);

  // Matrix of the quadratic form
  matrix<double> O_mat(J, J);

  int iter;
  for(iter = 0; iter < p.max_iter; ++iter) {
    if(p.verbosity >= 1)
      mpi_cout(comm) << "  Iteration " << (iter + 1) << " out of " << p.max_iter
                     << std::endl;

    // Update matrix of the quadratic form
    update_O_mat(A_block,
                 Ap_block,
                 App_block,
                 p.unity_sum_coeff,
                 Q_k,
                 D_k,
                 B_k,
                 p.default_model_weights,
                 index_map,
                 comm,
                 p.verbosity >= 2,
                 O_mat);

    // Minimize the functional
    c = inverse(O_mat) * b;

    //
    // Update regularization parameters
    //

    DD *= p.der_penalty_coeff;
    D_k *= p.der_penalty_coeff;
    B_k *= p.der_penalty_coeff;

    triqs::clef::placeholder<0> j_;
    triqs::clef::placeholder<1> k_;

    // Linear combination of particular solutions: Rank-local stage
    if(my_J != 0) {
      A_k(k_) << triqs::clef::sum(c(my_j_start + j_) * A_block(j_, k_),
                                  j_ = range(my_j_size));
      Ap_k(k_) << triqs::clef::sum(c(my_j_start + j_) * Ap_block(j_, k_),
                                   j_ = range(my_j_size));
      App_k(k_) << triqs::clef::sum(c(my_j_start + j_) * App_block(j_, k_),
                                    j_ = range(my_j_size));
    } else {
      A_k() = 0;
      Ap_k() = 0;
      App_k() = 0;
    }

    // Linear combination of particular solutions: Inter-rank stage
    A_k = mpi::all_reduce(A_k, comm);
    Ap_k = mpi::all_reduce(Ap_k, comm);
    App_k = mpi::all_reduce(App_k, comm);

    // Update Q
    for(int k : range(A_k.size())) {
      if(A_k(k) < 0)
        Q_k(k) = p.amp_penalty_max;
      else
        Q_k(k) /= p.amp_penalty_divisor;
    }

    // Update D
    for(int k : range(Ap_k.size())) {
      double val = std::abs(Ap_k(k));
      if(D_k(k) * val > DD) D_k(k) = DD / val;
    }

    // Update B
    for(int k : range(App_k.size())) {
      double val = std::abs(App_k(k));
      if(B_k(k) * val > DD) B_k(k) = DD / val;
    }

    double diff = std::abs(sum(c) - 1);
    if(p.verbosity >= 1) {
      double sum_abs = sum(abs(c));
      mpi_cout(comm) << "  End of iteration " << (iter + 1)
                     << ": |sum(c) - 1| = " << diff
                     << ", sum(abs(c)) = " << sum_abs << std::endl;
    }

    if(diff < convergence_tol) {
      if(p.verbosity >= 1) mpi_cout(comm) << "Convergence reached" << std::endl;
      break;
    }
  }

  if(p.verbosity >= 1 && iter == p.max_iter)
    mpi_cout(comm) << "Maximum number of iterations reached" << std::endl;

  return c;
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

  for(int n : range(data.size())) {
    auto& d = data[n];

    if(p.verbosity > 0)
      mpi_cout(comm) << "Using CC protocol to construct the final solution "
                        "for observable component ["
                     << n << "," << n << "]" << std::endl;

    // Select good solutions
    std::vector<int> good_solution_indices;
    good_solution_indices.reserve(d.particular_solutions.size());
    double chi_min = std::sqrt(d.objf_min);
    for(auto const& [s_index, s] : enumerate(d.particular_solutions)) {
      double chi = std::sqrt(s.second);
      if(chi / chi_min <= p.good_chi_rel && chi <= p.good_chi_abs)
        good_solution_indices.emplace_back(s_index);
    }

    if(p.verbosity > 1)
      mpi_cout(comm) << "Selected " << good_solution_indices.size()
                     << " good solutions" << std::endl;

    global_index_map index_map(comm, good_solution_indices.size());

    // Convergence tolerance for CC iterations
    using mesh_t = typename KernelType::mesh_type;
    double convergence_tol =
        min_element(abs(d.get_error_bars<mesh_t>() / d.get_rhs<mesh_t>()));
    if(p.verbosity > 1)
      mpi_cout(comm) << "Convergence tolerance is " << convergence_tol
                     << std::endl;

    // Compute coefficients c_j using CC optimization protocol
    auto c = cc_protocol(d.particular_solutions,
                         good_solution_indices,
                         e_k_list,
                         convergence_tol,
                         p,
                         index_map,
                         comm);

    if(p.verbosity >= 1) {
      mpi_cout(comm)
          << "Forming the resulting linear combination of particular solutions"
          << std::endl;
      if(p.verbosity >= 2)
        mpi_cout(comm) << "Coefficients of the linear combination: " << c
                       << std::endl;
    }

    configuration sol_sum(ci);
    int const j_range_start = index_map.range_start(comm.rank());
    for(auto const& [j, s_index] : enumerate(good_solution_indices)) {
      sol_sum += c(j_range_start + j) * d.particular_solutions[s_index].first;
    }
    d.final_solution = mpi::all_reduce(sol_sum, comm);
  }

  return compute_objf_final();
}

std::vector<double>
som_core::compute_final_solution_cc(final_solution_cc_parameters_t const& p) {
  p.validate();

#define IMPL_CASE(r, okmk)                                                     \
  case(kernel_id(BOOST_PP_SEQ_ELEM(0, okmk), BOOST_PP_SEQ_ELEM(1, okmk){})):   \
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
