/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016 by I. Krivenko
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

#include <vector>

#include <triqs/arrays/vector.hpp>
#include <triqs/gfs.hpp>
#include <triqs/utility/variant.hpp>
#include <triqs/mpi/base.hpp>
#include <triqs/statistics/histograms.hpp>

#include "run_parameters.hpp"
#include "configuration.hpp"
#include "kernels/base.hpp"

namespace som {

using namespace triqs::arrays;
using namespace triqs::gfs;
using triqs::utility::variant;
using triqs::statistics::histogram;

class som_core {

 // LHS cache index
 // Must be sure this object is destroyed after all other members
 cache_index ci;

 using input_data_r_t = std::vector<vector<double>>;
 using input_data_c_t = std::vector<vector<std::complex<double>>>;

 template<typename Mesh> using input_data_t = std14::conditional_t<
  std::is_same<Mesh,gf_mesh<imfreq>>::value,input_data_c_t,input_data_r_t>;

 // The right-hand side of the Fredholm integral equation
 // One vector per diagonal matrix element of the Green's function
 variant<input_data_r_t,input_data_c_t> rhs;
 // Error bars of the RHS, see eq. 30
 // One vector per diagonal matrix element of the Green's function
 variant<input_data_r_t,input_data_c_t> error_bars;

 // Resulting configurations
 std::vector<configuration> results;

 // Kind of the observable, GF/Susceptibility/Conductivity
 observable_kind kind;

 // Mesh of the input functions
 variant<gf_mesh<imtime>, gf_mesh<imfreq>, gf_mesh<legendre>> mesh;

 // Norm of the solution to be found
 double norm;

 // Parameters of the last call to run()
 run_parameters_t params;

 // Status of the run upon exit: 0 for clean termination, > 0 otherwise.
 int run_status = 0;

 // MPI communicator
 triqs::mpi::communicator comm;

 // Objective function histograms
 std::vector<histogram> histograms;

 // Fill rhs and error_bars
 template<typename... GfOpts>
 void set_input_data(gf_const_view<GfOpts...> g, gf_const_view<GfOpts...> S);

 // Run the main part of the algorithm
 template<typename KernelType> void run_impl();

 // Adjust the number of global updates (F)
 template<typename KernelType> int adjust_f(
  KernelType const& kern,
  typename KernelType::result_type rhs_,
  typename KernelType::result_type error_bars_,
  std::function<bool()> const& stop_callback);

 // Accumulate solutions
 template<typename KernelType> configuration accumulate(
  KernelType const& kern,
  typename KernelType::result_type rhs_,
  typename KernelType::result_type error_bars_,
  histogram & hist,
  std::function<bool()> const& stop_callback,
  int F);

public:

 /// Construct on imaginary-time quantities
 som_core(gf_const_view<imtime> g_tau, gf_const_view<imtime> S_tau,
          observable_kind kind = FermionGf, double norm = 1.0);
 /// Construct on imaginary-frequency quantities
 som_core(gf_const_view<imfreq> g_iw, gf_const_view<imfreq> S_iw,
          observable_kind kind = FermionGf, double norm = 1.0);
 /// Construct on quantities in Legendre polynomial basis
 som_core(gf_const_view<legendre> g_l, gf_const_view<legendre> S_l,
          observable_kind kind = FermionGf, double norm = 1.0);

 /// Destructor
 ~som_core();

 TRIQS_WRAP_ARG_AS_DICT // Wrap the parameters as a dictionary in python with c++2py
 void run(run_parameters_t const& p);

 /// Set of parameters used in the last call to run()
 run_parameters_t get_last_run_parameters() const { return params; }

 /// Status of the run on exit
 int get_run_status() const { return run_status; }

 /// Fill a Green's function using the calculated spectra
 template<typename MeshType>
 friend void triqs_gf_view_assign_delegation(gf_view<MeshType> g, som_core const& cont);

 /// Accumulated objective function histograms
 std::vector<histogram> const& get_histograms() const { return histograms; }

};

}
