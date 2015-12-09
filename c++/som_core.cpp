#include "som_core.hpp"
#include "global_counter.hpp"

namespace som {

template<> int som_core::get_mesh_id<imtime>() { return 0; }
template<> int som_core::get_mesh_id<imfreq>() { return 1; }
template<> int som_core::get_mesh_id<legendre>() { return 2; }

template<typename... GfOpts>
void check_input_gf(gf_const_view<GfOpts...> g, gf_const_view<GfOpts...> S) {
 if(g.mesh() != S.mesh() || get_target_shape(g) != get_target_shape(S))
  TRIQS_RUNTIME_ERROR << "som_core: input quantity and the error-bar function S(\\tau) must have equivalent structure";

 auto shape = get_target_shape(g);
 if(shape[0] != shape[1]) TRIQS_RUNTIME_ERROR << "som_core: matrix-valued input quantities must be square";
}

template<typename... GfOpts>
void som_core::set_input_data(gf_const_view<GfOpts...> g, gf_const_view<GfOpts...> S) {
 using mesh_t = typename gf_const_view<GfOpts...>::mesh_t;
 using input_data_t = std14::conditional_t<std::is_same<mesh_t,gf_mesh<imfreq>>::value,input_data_c_t,input_data_r_t>;

 auto & rhs_ = (input_data_t&)(rhs);
 auto & error_bars_ = (input_data_t&)error_bars;
 auto gf_dim = get_target_shape(g)[0];
 for(int i = 0; i < gf_dim; ++i) {
  rhs_.emplace_back(g.mesh().size());
  rhs_.back()() = g.data()(range(),i,i);

  error_bars_.emplace_back(S.mesh().size());
  error_bars_.back()() = S.data()(range(),i,i);
 }

 results.resize(gf_dim);
}

////////////////////////////////
// Constructors: fermionic GF //
////////////////////////////////

// Construct on fermionic G(\tau)
som_core::som_core(gf_const_view<imtime> g_tau, gf_const_view<imtime> S_tau) :
kind(FermionicGf), mesh_id(get_mesh_id<imtime>()),
rhs(input_data_r_t()), error_bars(input_data_r_t()) {
 if(g_tau.domain().statistic != Fermion)
  TRIQS_RUNTIME_ERROR << "som_core: only fermionic Green's functions are supported";
 check_input_gf(g_tau,S_tau);
 if(!is_gf_real(g_tau) || !is_gf_real(S_tau))
  TRIQS_RUNTIME_ERROR << "som_core: imaginary time Green's functions must be real";
 gf<imtime, matrix_real_valued> g_tau_real = real(g_tau), S_tau_real = real(S_tau);
 set_input_data(make_const_view(g_tau_real),make_const_view(S_tau_real));
}

// Construct on fermionic G(i\omega)
som_core::som_core(gf_const_view<imfreq> g_iw, gf_const_view<imfreq> S_iw) :
kind(FermionicGf), mesh_id(get_mesh_id<imfreq>()),
rhs(input_data_c_t()), error_bars(input_data_c_t()) {
 if(g_iw.domain().statistic != Fermion)
  TRIQS_RUNTIME_ERROR << "som_core: only fermionic Green's functions are supported";
 check_input_gf(g_iw,S_iw);
 set_input_data(g_iw,S_iw);
}

// Construct on fermionic G_l
som_core::som_core(gf_const_view<legendre> g_l, gf_const_view<legendre> S_l) :
kind(FermionicGf), mesh_id(get_mesh_id<legendre>()),
rhs(input_data_r_t()), error_bars(input_data_r_t()) {
 if(g_l.domain().statistic != Fermion)
  TRIQS_RUNTIME_ERROR << "som_core: only fermionic Green's functions are supported";
 check_input_gf(g_l,S_l);
 if(!is_gf_real(g_l) || !is_gf_real(S_l))
  TRIQS_RUNTIME_ERROR << "som_core: Legendre Green's functions must be real";
 gf<legendre, matrix_real_valued> g_l_real = real(g_l), S_l_real = real(S_l);
 set_input_data(make_const_view(g_l_real),make_const_view(S_l_real));
}

void som_core::run(run_parameters_t const& p) {
 // FIXME
 for(auto & r :results) r.insert({2.0,2.0,1.5});
}

void triqs_gf_view_assign_delegation(gf_view<refreq> g_w, som_core const& cont) {
 auto shape = get_target_shape(g_w);
 auto gf_dim = cont.results.size();

 if(shape[0] != gf_dim || shape[1] != gf_dim)
  TRIQS_RUNTIME_ERROR << "som_core: expected a real-frequency Green's function with matrix dimensions "
                      << gf_dim << "x" << gf_dim << " in assignment";

 g_w() = 0;
 for(int i = 0; i < gf_dim; ++i) {
  auto const& conf = cont.results[i];
  for(auto const& rect : conf) {
   for(auto e : g_w.mesh()) g_w.data()(e.index(),i,i) += rect.hilbert_transform(double(e));
   auto & tail = g_w.singularity();
   tail.data()(range(),i,i) += rect.tail_coefficients(tail.order_min(),tail.order_max());
  }
 }
}

}
