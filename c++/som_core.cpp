#include "som_core.hpp"
#include "global_counter.hpp"

namespace som {

som_core::som_core(gf_const_view<imtime> g_tau, gf_const_view<imtime> S) {
 if(g_tau.domain().statistic != Fermion)
  TRIQS_RUNTIME_ERROR << "som_core: only fermionic Green's functions are supported";
 if(g_tau.mesh() != S.mesh() || get_target_shape(g_tau) != get_target_shape(S))
  TRIQS_RUNTIME_ERROR << "som_core: G(\\tau) and the error-bar function S(\\tau) must have equivalent structure";
}

void som_core::run(run_parameters_t const& p) {
}

}
