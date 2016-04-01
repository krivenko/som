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
 cache_index ci;
 kernel<FermionGf,imtime> kern(mesh);
 obj_function of(kern, g_tau, s_tau);

 auto params = run_parameters_t({-3.0,3.0});
 params.n_elementary_updates = 100;
 params.min_rect_width = 0.001;
 params.min_rect_weight = 0.001;

 solution_worker<kernel<FermionGf,imtime>> worker(of, 1.0, ci, params, 10);

 auto solution = worker(10);

 configuration solution_ref(ci);
 triqs::h5::file arch("solution_worker.ref.h5", H5F_ACC_RDONLY);
 h5_read(arch, "RandomConfig_output", solution_ref, ci);

 EXPECT_EQ(solution_ref, solution);
}

TEST(solution_worker,StartConfig) {
 cache_index ci;
 kernel<FermionGf,imtime> kern(mesh);
 obj_function of(kern, g_tau, s_tau);

 auto params = run_parameters_t({-3.0,3.0});
 params.n_elementary_updates = 100;
 params.min_rect_width = 0.001;
 params.min_rect_weight = 0.001;

 solution_worker<kernel<FermionGf,imtime>> worker(of, 1.0, ci, params, 10);

 triqs::h5::file arch("solution_worker.ref.h5", H5F_ACC_RDONLY);
 configuration init_config(ci);
 h5_read(arch, "StartConfig_input", init_config, ci);

 auto solution = worker(init_config);

 configuration solution_ref(ci);
 h5_read(arch, "StartConfig_output", solution_ref, ci);

 EXPECT_EQ(solution_ref, solution);
}

MAKE_MAIN;
