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
#pragma once

#include <vector>

#include <triqs/arrays/vector.hpp>
#include <triqs/gfs.hpp>
#include <triqs/utility/variant.hpp>

#include "run_parameters.hpp"
#include "configuration.hpp"
#include "kernels/fermiongf_imtime.hpp"
#include "objective_function.hpp"

namespace som {

using namespace triqs::arrays;
using namespace triqs::gfs;
using triqs::utility::variant;

class som_core {

 using input_data_r_t = std::vector<vector<double>>;
 using input_data_c_t = std::vector<vector<std::complex<double>>>;

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

 // Mesh type, imtime/imfreq/legendre
 int mesh_id;
 template<typename MT> static int get_mesh_id();

 // Fill rhs and error_bars
 template<typename... GfOpts>
 void set_input_data(gf_const_view<GfOpts...> g, gf_const_view<GfOpts...> S);

public:

 som_core(gf_const_view<imtime> g_tau, gf_const_view<imtime> S_tau);
 som_core(gf_const_view<imfreq> g_iw, gf_const_view<imfreq> S_iw);
 som_core(gf_const_view<legendre> g_l, gf_const_view<legendre> S_l);

 TRIQS_WRAP_ARG_AS_DICT // Wrap the parameters as a dictionary in python with the clang tool
 void run(run_parameters_t const& p);

 gf_view<refreq> operator()(gf_view<refreq> g_w) const;

 // Fill gf<refreq> with obtained results
 friend void triqs_gf_view_assign_delegation(gf_view<refreq> g_w, som_core const& cont) { cont(g_w); }

};

}
