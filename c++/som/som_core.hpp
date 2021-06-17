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
#pragma once

#include <complex>
#include <optional>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <mpi/mpi.hpp>
#include <triqs/arrays/vector.hpp>
#include <triqs/gfs.hpp>
#include <triqs/statistics/histograms.hpp>

#include "configuration.hpp"
#include "kernels/observables.hpp"
#include "run_parameters.hpp"

namespace som {

class som_core {

  // LHS cache index
  // Must be sure this object is destroyed after all other members
  cache_index ci;

  // Kind of the observable, GF/Susceptibility/Conductivity
  observable_kind kind;

  // Mesh of the input functions
  std::variant<triqs::gfs::gf_mesh<triqs::gfs::imtime>,
               triqs::gfs::gf_mesh<triqs::gfs::imfreq>,
               triqs::gfs::gf_mesh<triqs::gfs::legendre>>
      mesh;

  using histogram_t = triqs::statistics::histogram;

  // Input/output data (one structure per diagonal matrix element of g)
  struct data_t {

    using input_data_r_t = triqs::arrays::array<double, 1>;
    using input_data_c_t = triqs::arrays::array<std::complex<double>, 1>;

    template <typename Mesh>
    using input_data_t = std::conditional_t<
        std::is_same<Mesh, triqs::gfs::gf_mesh<triqs::gfs::imfreq>>::value,
        input_data_c_t, input_data_r_t>;

    // The right-hand side of the Fredholm integral equation
    std::variant<input_data_r_t, input_data_c_t> rhs;
    // Error bars of the RHS
    std::variant<input_data_r_t, input_data_c_t> error_bars;

    // Norm of the solutions to be found
    double norm = 1.0;

    // Accumulated pairs (basis solution, objective function)
    std::vector<std::pair<configuration, double>> basis_solutions;

    // Final solution
    configuration final_solution;

    // Objective function histogram
    std::optional<histogram_t> histogram;

    // Constructor
    template <typename Mesh> data_t(Mesh const& mesh, cache_index & ci);

    // Initialize 'rhs', 'error_bars' and 'norm'
    template <typename... GfOpts>
    void init_input(int i,
                    triqs::gfs::gf_const_view<GfOpts...> g,
                    triqs::gfs::gf_const_view<GfOpts...> S,
                    double norm_
                    );

    // Typed access to 'rhs'
    template <typename Mesh> input_data_t<Mesh> const& get_rhs() const {
      return std::get<input_data_t<Mesh>>(rhs);
    }
    // Typed access to 'error_bars'
    template <typename Mesh> input_data_t<Mesh> const& get_error_bars() const {
      return std::get<input_data_t<Mesh>>(error_bars);
    }
  };
  std::vector<data_t> data;

  // Parameters of the last call to run()
  run_parameters_t params;

  // Status of the run upon exit: 0 for clean termination, > 0 otherwise.
  int run_status = 0;

  // MPI communicator
  mpi::communicator comm;

  // Run the main part of the algorithm
  template <typename KernelType> void run_impl();

  // Adjust the number of global updates (F)
  template <typename KernelType>
  int adjust_f(KernelType const& kern,
               data_t const& data,
               std::function<bool()> const& stop_callback);

  // Accumulate solutions
  template <typename KernelType>
  configuration accumulate(KernelType const& kern,
                           data_t & data,
                           std::function<bool()> const& stop_callback, int F);

public:
  /// Construct on imaginary-time quantities
  som_core(triqs::gfs::gf_const_view<triqs::gfs::imtime> g_tau,
           triqs::gfs::gf_const_view<triqs::gfs::imtime> S_tau,
           observable_kind kind = FermionGf,
           triqs::arrays::vector<double> const& norms = {});
  /// Construct on imaginary-frequency quantities
  som_core(triqs::gfs::gf_const_view<triqs::gfs::imfreq> g_iw,
           triqs::gfs::gf_const_view<triqs::gfs::imfreq> S_iw,
           observable_kind kind = FermionGf,
           triqs::arrays::vector<double> const& norms = {});
  /// Construct on quantities in Legendre polynomial basis
  som_core(triqs::gfs::gf_const_view<triqs::gfs::legendre> g_l,
           triqs::gfs::gf_const_view<triqs::gfs::legendre> S_l,
           observable_kind kind = FermionGf,
           triqs::arrays::vector<double> const& norms = {});

  // Wrap the parameters as a dictionary in python with c++2py
  TRIQS_WRAP_ARG_AS_DICT void run(run_parameters_t const& p);

  /// Set of parameters used in the last call to run()
  run_parameters_t get_last_run_parameters() const { return params; }

  /// Status of the run on exit
  int get_run_status() const { return run_status; }

  /// Fill a Green's function using the calculated spectra
  template <typename MeshType>
  friend void triqs_gf_view_assign_delegation(triqs::gfs::gf_view<MeshType> g,
                                              som_core const& cont);

  /// Compute GF tail coefficients using the calculated spectra
  [[nodiscard]] triqs::arrays::array<std::complex<double>, 3>
  compute_tail(int max_order) const;

  /// Matrix dimension of the observable
  [[nodiscard]] int get_dim() const { return data.size(); }

  /// Final solution for the i-th diagonal matrix element
  [[nodiscard]] configuration const& get_solution(int i) const;

  /// Final solutions
  [[nodiscard]] std::vector<configuration> get_solutions() const;

  /// Accumulated objective function histogram for the i-th diagonal matrix
  /// element of the observable
  [[nodiscard]] std::optional<histogram_t> const& get_histogram(int i) const;

  /// Accumulated objective function histogram for all diagonal matrix
  /// elements of the observable
  [[nodiscard]] std::optional<std::vector<histogram_t>> get_histograms() const;

  /// Discard all accumulated basis solutions, histograms and final solutions
  void clear();
};

} // namespace som
