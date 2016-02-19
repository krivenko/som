#include <triqs/test_tools/arrays.hpp>
#include <triqs/h5.hpp>

#include "kernels/fermiongf_imtime.hpp"
#include "solution_worker.hpp"

using namespace som;
using namespace triqs::arrays;

double beta;
gf_mesh<imtime> mesh;
array<double,1> g_tau, s_tau;

struct Env : public ::testing::Environment {
 virtual void SetUp() {
  triqs::h5::file arch("solution_worker.ref.h5", H5F_ACC_RDONLY);

  h5_read(arch, "beta", beta);
  h5_read(arch, "g_tau", g_tau);
  h5_read(arch, "s_tau", s_tau);

  mesh = {beta, Fermion, static_cast<int>(first_dim(g_tau))};
 }
};
::testing::Environment* const env = ::testing::AddGlobalTestEnvironment(new Env);

using obj_function = objective_function<kernel<FermionGf,imtime>>;

TEST(solution_worker,RandomConfig) {
 obj_function of(mesh, g_tau, s_tau);

 auto params = run_parameters_t({-3.0,3.0});
 params.n_elementary_updates = 100;
 params.min_rect_width = 0.001;
 params.min_rect_weight = 0.001;

 solution_worker<kernel<FermionGf,imtime>> worker(of, 1.0, params, 10);

 worker(10);

 // TODO: compare to the reference solution
}

TEST(solution_worker,StartConfig) {
 // TODO
}

MAKE_MAIN;
