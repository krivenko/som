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
#pragma once

#include <cmath>
#include <complex>
#include <limits>
#include <optional>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <mpi/mpi.hpp>
#include <nda/nda.hpp>
#include <triqs/gfs.hpp>
#include <triqs/stat/histograms.hpp>

#include "parameters.hpp"
#include <som/configuration.hpp>
#include <som/kernels/mesh_traits.hpp>
#include <som/kernels/observables.hpp>

namespace testing {
class spectral_stats_test;
}

namespace som {

class som_core {

  // For unit testing
  friend class testing::spectral_stats_test;
  som_core() = default;

  // LHS cache index
  // Must be sure this object is destroyed after all other members
  cache_index ci;

  // Kind of the observable, GF/Susceptibility/Conductivity
  observable_kind kind = {};

  // Mesh of the input functions
  mesh_variant_t mesh;

  using histogram_t = triqs::stat::histogram;

  // Input/output data (one structure per diagonal matrix element of g)
  struct data_t {

    using input_data_r_t = nda::array<double, 1>;
    using input_data_c_t = nda::array<std::complex<double>, 1>;

    template <typename Mesh>
    using input_data_t =
        std::conditional_t<std::is_same<Mesh, triqs::mesh::imfreq>::value,
                           input_data_c_t,
                           input_data_r_t>;

    // The right-hand side of the Fredholm integral equation
    std::variant<input_data_r_t, input_data_c_t> rhs;
    // Error bars of the RHS
    std::variant<input_data_r_t, input_data_c_t> error_bars;

    // Norm of the solutions to be found
    double norm = 1.0;

    // Accumulated pairs (particular solution, objective function).
    // For the sake of memory efficiency, this list of solutions is rank-local.
    std::vector<std::pair<configuration, double>> particular_solutions;

    // Minimum of :math:`\chi^2` over all particular solutions
    double objf_min = HUGE_VAL;

    // Final solution
    configuration final_solution;

    // :math:`\chi^2` of the final solution
    double objf_final = std::numeric_limits<double>::quiet_NaN();

    // Objective function histogram
    std::optional<histogram_t> histogram;

    // Constructor
    template <typename Mesh> data_t(Mesh const& mesh, cache_index& ci);

    // Initialize 'rhs', 'error_bars' and 'norm'
    template <typename... GfOpts>
    void init_input(long i,
                    triqs::gfs::gf_const_view<GfOpts...> g,
                    triqs::gfs::gf_const_view<GfOpts...> S,
                    double norm_);

    // Typed access to 'rhs'
    template <typename Mesh> input_data_t<Mesh> const& get_rhs() const {
      return std::get<input_data_t<Mesh>>(rhs);
    }
    // Typed access to 'error_bars'
    template <typename Mesh> input_data_t<Mesh> const& get_error_bars() const {
      return std::get<input_data_t<Mesh>>(error_bars);
    }

    template <typename KernelType>
    void compute_objf_final(KernelType const& kernel);
  };
  std::vector<data_t> data;

  // Parameters of the last call to accumulate()
  accumulate_parameters_t params;

  // Status of accumulate() upon exit: 0 for clean termination, > 0 otherwise.
  int accumulate_status = 0;

  // MPI communicator
  mpi::communicator comm;

  // Implementation details of adjust_f()
  template <typename KernelType>
  int adjust_f_impl(adjust_f_parameters_t const& params);

  // Implementation details of accumulate()
  template <typename KernelType> void accumulate_impl();

  // Compute objective function :math:`\chi^2` of final solutions
  std::vector<double> compute_objf_final();
  template <typename KernelType> std::vector<double> compute_objf_final_impl();

  // Implementation details of compute_final_solution_cc_impl()
  template <typename KernelType>
  std::vector<double>
  compute_final_solution_cc_impl(final_solution_cc_parameters_t const& params);

public:
  /// Construct on imaginary-time quantities
  som_core(triqs::gfs::gf_const_view<triqs::mesh::imtime> g_tau,
           triqs::gfs::gf_const_view<triqs::mesh::imtime> S_tau,
           observable_kind kind = FermionGf,
           nda::vector<double> const& norms = {});
  /// Construct on imaginary-frequency quantities
  som_core(triqs::gfs::gf_const_view<triqs::mesh::imfreq> g_iw,
           triqs::gfs::gf_const_view<triqs::mesh::imfreq> S_iw,
           observable_kind kind = FermionGf,
           nda::vector<double> const& norms = {});
  /// Construct on quantities in Legendre polynomial basis
  som_core(triqs::gfs::gf_const_view<triqs::mesh::legendre> g_l,
           triqs::gfs::gf_const_view<triqs::mesh::legendre> S_l,
           observable_kind kind = FermionGf,
           nda::vector<double> const& norms = {});

  // Automatically adjust the number of global updates (F)
  TRIQS_WRAP_ARG_AS_DICT int adjust_f(adjust_f_parameters_t const& p);

  /// Accumulate particular solutions
  TRIQS_WRAP_ARG_AS_DICT void accumulate(accumulate_parameters_t const& p);

  /// Select particular solutions according to the standard SOM criterion
  /// and compute the final solution. Selected solutions must fulfill
  /// :math:`\chi/\chi_{min} \leq good_chi_rel` *and*
  /// :math:`\chi \leq good_chi_abs`.
  std::vector<double> compute_final_solution(double good_chi_rel = 2.0,
                                             double good_chi_abs = HUGE_VAL,
                                             int verbosity = 0);

  /// Compute the final solution using the SOCC protocol
  TRIQS_WRAP_ARG_AS_DICT
  std::vector<double>
  compute_final_solution_cc(final_solution_cc_parameters_t const& p);

  /// Set of parameters used in the last call to accumulate()
  [[nodiscard]] accumulate_parameters_t get_last_accumulate_parameters() const {
    return params;
  }

  /// Status of the accumulate() on exit
  [[nodiscard]] int get_accumulate_status() const { return accumulate_status; }

  /// Accumulate particular solutions and compute the final solution using
  /// the standard SOM criterion
  TRIQS_WRAP_ARG_AS_DICT void run(accumulate_parameters_t const& p);

  /// Kind of the observable
  [[nodiscard]] observable_kind get_observable_kind() const { return kind; }

  /// Matrix dimension of the observable
  [[nodiscard]] long get_dim() const { return long(data.size()); }

  /// MPI communicator
  [[nodiscard]] mpi::communicator const& get_comm() const { return comm; }

  /// Accumulated particular solutions for the i-th diagonal matrix element.
  /// The returned list of solutions is rank-local.
  [[nodiscard]] std::vector<std::pair<configuration, double>> const&
  get_particular_solutions(long i) const;

  /// Final solution for the i-th diagonal matrix element
  [[nodiscard]] configuration const& get_solution(long i) const;

  /// Final solutions
  [[nodiscard]] std::vector<configuration> get_solutions() const;

  /// Value of the objective function :math:`\chi^2` of the final solution
  /// for the i-th diagonal matrix element.
  [[nodiscard]] double get_objf(long i) const;

  /// Values of the objective function :math:`\chi^2` of the final solutions
  /// (one value per a diagonal matrix element of the observable).
  [[nodiscard]] std::vector<double> get_objf() const;

  /// Accumulated objective function histogram for the i-th diagonal matrix
  /// element of the observable
  [[nodiscard]] std::optional<histogram_t> const& get_histogram(long i) const;

  /// Accumulated objective function histogram for all diagonal matrix
  /// elements of the observable
  [[nodiscard]] std::optional<std::vector<histogram_t>> get_histograms() const;

  /// Minimum of the objective function over all accumulated particular
  /// solutions (one value per a diagonal matrix element of the observable)
  [[nodiscard]] std::vector<double> get_objf_min() const;

  /// Discard all accumulated particular solutions, histograms
  /// and final solutions
  void clear();
};

/// Fill a real-frequency observable from a computed SOM solution
void fill_refreq(triqs::gfs::gf_view<triqs::mesh::refreq> g_w,
                 som_core const& cont,
                 bool with_binning = false);
/// Fill a real-frequency observable from a list of solutions
/// (one solution per a diagonal matrix element of the observable)
void fill_refreq(triqs::gfs::gf_view<triqs::mesh::refreq> g_w,
                 observable_kind kind,
                 std::vector<configuration> const& solutions,
                 bool with_binning = false);

/// Compute tail coefficients from a computed SOM solution
[[nodiscard]] nda::array<std::complex<double>, 3>
compute_tail(int max_order, som_core const& cont);
/// Compute tail coefficients from a list of solutions
/// (one solution per a diagonal matrix element of the observable)
[[nodiscard]] nda::array<std::complex<double>, 3>
compute_tail(int max_order,
             observable_kind kind,
             std::vector<configuration> const& solutions);

/// Reconstruct an observable in the imaginary-time representation
/// from a computed SOM solution
void reconstruct(triqs::gfs::gf_view<triqs::mesh::imtime> g,
                 som_core const& cont);
/// Reconstruct an observable in the imaginary-frequency representation
/// from a computed SOM solution
void reconstruct(triqs::gfs::gf_view<triqs::mesh::imfreq> g,
                 som_core const& cont);
/// Reconstruct an observable in the Legendre polynomial basis
/// from a computed SOM solution
void reconstruct(triqs::gfs::gf_view<triqs::mesh::legendre> g,
                 som_core const& cont);

/// Reconstruct an observable in the imaginary-time representation from a list
/// of solutions (one solution per a diagonal matrix element of the observable)
void reconstruct(triqs::gfs::gf_view<triqs::mesh::imtime> g,
                 observable_kind kind,
                 std::vector<configuration> const& solutions);
/// Reconstruct an observable in the imaginary-frequency representation from
/// a list of solutions (one solution per a diagonal matrix element of the
/// observable)
void reconstruct(triqs::gfs::gf_view<triqs::mesh::imfreq> g,
                 observable_kind kind,
                 std::vector<configuration> const& solutions);
/// Reconstruct an observable in the Legendre polynomial basis from a list of
/// solutions (one solution per a diagonal matrix element of the observable)
void reconstruct(triqs::gfs::gf_view<triqs::mesh::legendre> g,
                 observable_kind kind,
                 std::vector<configuration> const& solutions);

} // namespace som
