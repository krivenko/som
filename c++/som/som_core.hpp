/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
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
#include <vector>

#include <mpi/mpi.hpp>
#include <triqs/arrays/vector.hpp>
#include <triqs/gfs.hpp>
#include <triqs/statistics/histograms.hpp>
#include <triqs/utility/variant.hpp>

#include "configuration.hpp"
#include "kernels/observables.hpp"
#include "run_parameters.hpp"

namespace som {

class som_core {

  // LHS cache index
  // Must be sure this object is destroyed after all other members
  cache_index ci;

  using input_data_r_t = std::vector<triqs::arrays::vector<double>>;
  using input_data_c_t =
      std::vector<triqs::arrays::vector<std::complex<double>>>;

  template <typename Mesh>
  using input_data_t = std::conditional_t<
      std::is_same<Mesh, triqs::gfs::gf_mesh<triqs::gfs::imfreq>>::value,
      input_data_c_t, input_data_r_t>;

  // Kind of the observable, GF/Susceptibility/Conductivity
  observable_kind kind;

  // Mesh of the input functions
  std::variant<triqs::gfs::gf_mesh<triqs::gfs::imtime>,
               triqs::gfs::gf_mesh<triqs::gfs::imfreq>,
               triqs::gfs::gf_mesh<triqs::gfs::legendre>>
      mesh;

  // The right-hand side of the Fredholm integral equation
  // One vector per diagonal matrix element of the Green's function
  std::variant<input_data_r_t, input_data_c_t> rhs;
  // Error bars of the RHS, see eq. 30
  // One vector per diagonal matrix element of the Green's function
  std::variant<input_data_r_t, input_data_c_t> error_bars;

  // Norms of the solutions to be found (one number per diagonal matrix element
  // of g)
  triqs::arrays::vector<double> norms;

  // Resulting configurations
  std::vector<configuration> results;

  // Parameters of the last call to run()
  run_parameters_t params;

  // Status of the run upon exit: 0 for clean termination, > 0 otherwise.
  int run_status = 0;

  // MPI communicator
  mpi::communicator comm;

  // Objective function histograms
  std::vector<triqs::statistics::histogram> histograms;

  // Fill rhs and error_bars
  template <typename... GfOpts>
  void set_input_data(triqs::gfs::gf_const_view<GfOpts...> g,
                      triqs::gfs::gf_const_view<GfOpts...> S);

  // Run the main part of the algorithm
  template <typename KernelType> void run_impl();

  // Adjust the number of global updates (F)
  template <typename KernelType>
  int adjust_f(KernelType const& kern, typename KernelType::result_type rhs_,
               typename KernelType::result_type error_bars_, double norm,
               std::function<bool()> const& stop_callback);

  // Accumulate solutions
  template <typename KernelType>
  configuration accumulate(KernelType const& kern,
                           typename KernelType::result_type rhs_,
                           typename KernelType::result_type error_bars_,
                           double norm, triqs::statistics::histogram& hist,
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

  /// Accumulated solutions
  [[nodiscard]] std::vector<configuration> const& get_solutions() const {
    return results;
  }

  /// Accumulated objective function histograms
  [[nodiscard]] std::vector<triqs::statistics::histogram> const&
  get_histograms() const {
    return histograms;
  }
};

} // namespace som
