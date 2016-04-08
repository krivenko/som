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
#include "som_core.hpp"
#include "global_counter.hpp"

namespace som {

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
 using input_data_t = std14::conditional_t<std::is_same<mesh_t,gf_mesh<imfreq>>::value,
                                           input_data_c_t,input_data_r_t>;

 auto & rhs_ = (input_data_t&)(rhs);
 auto & error_bars_ = (input_data_t&)error_bars;
 auto gf_dim = get_target_shape(g)[0];

 results.reserve(gf_dim);
 for(int i = 0; i < gf_dim; ++i) {
  rhs_.emplace_back(g.data()(range(),i,i));
  error_bars_.emplace_back(S.data()(range(),i,i));
  results.emplace_back(ci);
 }
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

 auto & lp = last_run_parameters; lp = p;

 if((kind == Susceptibility || kind == Conductivity) && lp.energy_window.first < 0) {
  last_run_parameters.energy_window.first = 0;
  if(lp.verbosity > 0) std::cout << "WARNING: left boundary of the energy window is reset to 0";
 }

 if(lp.energy_window.first >= lp.energy_window.second)
  fatal_error("wrong energy window [" + to_string(lp.energy_window.first) +
              ";" + to_string(lp.energy_window.second) + "]");

 double window_width = lp.energy_window.second - lp.energy_window.first;
 if(lp.min_rect_width <= 0 || lp.min_rect_width >= window_width)
  fatal_error("min_rect_width must be in [0;" + to_string(window_width) + "]");

 if(lp.min_rect_weight <= 0 || lp.min_rect_weight >= norm)
  fatal_error("min_rect_weight must be in [0;" + to_string(norm) + "]");

#define EI(ok, mk) int(ok) + 3 * mk

 switch(EI(kind, mesh.index())) {
  case EI(FermionGf,0): engine<kernel<FermionGf,imtime>>(); break;
  case EI(FermionGf,1): engine<kernel<FermionGf,imfreq>>(); break;
  //case EI(FermionGf,2): engine<kernel<FermionGf,legendre>>(); break;
  // TODO
 }
#undef EI
}

gf_view<refreq> som_core::operator()(gf_view<refreq> g_w) const {
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
 return g_w;
}

template<typename KernelType> void som_core::engine() {

 typename KernelType::mesh_type const& m = mesh;

 if(last_run_parameters.verbosity > 0) {
  std::cout << KernelType::name() << std::endl;
  std::cout << "Mesh: " << m << std::endl;
  std::cout << "Constructing integral kernel... ";
 }
 KernelType kernel(m);
 if(last_run_parameters.verbosity > 0) std::cout << "done" << std::endl;

 // TODO
}

}
