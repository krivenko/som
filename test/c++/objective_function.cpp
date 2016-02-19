#include <vector>

#include <triqs/test_tools/arrays.hpp>
#include <triqs/gfs.hpp>

#include "kernels/fermiongf_imtime.hpp"
#include "objective_function.hpp"

using namespace som;
using namespace triqs::arrays;
using namespace triqs::gfs;

const double beta = 2;
std::vector<rectangle> rects {{-2,2.6, 0.3},
                              {1.3,2.6,0.6},
                              {-0.5,2.6,0.5},
                              {2,2.6,0.7}};

gf_mesh<imtime> mesh(beta,Fermion,11);

array<double,1> GF() {
 array<double,1> res(mesh.size());
 for(auto pt : mesh) {
  int i = pt.index();
  double tau = double(pt);
  res(i) = -0.5*(std::exp(-tau*1.3)/(1 + std::exp(-beta*1.3)) +
                 std::exp(tau*0.7)/(1 + std::exp(beta*0.7)));
 }
 return res;
}

array<double,1> S() {
 array<double,1> res(mesh.size());
 for(auto pt : mesh) {
  int i = pt.index();
  double tau = double(pt);
  res(i) = 0.05*std::exp(-std::abs(tau - 0.5*beta));
 }
 return res;
}

using obj_function = objective_function<kernel<FermionGf,imtime>>;

TEST(objective_function,Change) {
 configuration conf {rects[0],rects[2]};
 auto gf = GF();
 auto s = S();

 obj_function of(mesh, gf, s);
 auto const& kern = of.get_kernel();

 EXPECT_NEAR(172.784,of(conf),1e-3);

 config_update cu;
 cu.change_rectangle(1, rects[3]);
 EXPECT_NEAR(167.785,of(conf, cu),1e-3);

 cu.reset();
 EXPECT_NEAR(172.784,of(conf, cu),1e-3);
 cu.change_rectangle(0,rects[1]);
 EXPECT_NEAR(301.098,of(conf, cu),1e-3);

 kern.update_config_cache(conf, cu);
 apply(conf, cu);
 EXPECT_NEAR(301.098,of(conf),1e-3);
}

TEST(objective_function,Add) {
 configuration conf {rects[0],rects[1]};
 auto gf = GF();
 auto s = S();

 obj_function of(mesh, gf, s);
 auto const& kern = of.get_kernel();

 EXPECT_NEAR(173.149,of(conf),1e-3);

 config_update cu;
 cu.add_rectangle(rects[2]);
 EXPECT_NEAR(400.832,of(conf, cu),1e-3);

 cu.reset();
 EXPECT_NEAR(173.149,of(conf, cu),1e-3);
 cu.add_rectangle(rects[3]);
 EXPECT_NEAR(405.862,of(conf, cu),1e-3);

 kern.update_config_cache(conf, cu);
 apply(conf, cu);
 EXPECT_NEAR(405.862,of(conf),1e-3);
}

TEST(objective_function,Remove) {
 configuration conf {rects[0],rects[1],rects[2]};
 auto gf = GF();
 auto s = S();

 obj_function of(mesh, gf, s);
 auto const& kern = of.get_kernel();

 EXPECT_NEAR(400.832,of(conf),1e-3);

 config_update cu;
 cu.remove_rectangle(1);
 EXPECT_NEAR(172.784,of(conf, cu),1e-3);

 cu.reset();
 EXPECT_NEAR(400.832,of(conf, cu),1e-3);
 cu.remove_rectangle(2);
 EXPECT_NEAR(173.149,of(conf, cu),1e-3);

 kern.update_config_cache(conf, cu);
 apply(conf, cu);
 EXPECT_NEAR(173.149,of(conf),1e-3);
}

/*
TEST(objective_function,Multiple) {
 configuration conf {rects[0],rects[1]};
 auto gf = GF();
 auto s = S();

 obj_function of(mesh, gf, s);
 auto const& kern = of.get_kernel();

 EXPECT_NEAR(173.149,of(conf),1e-3);

 config_update cu;
 cu.add_rectangle(rects[3]);
 cu.change_rectangle(1,rects[2]);
 EXPECT_NEAR(395.468,of(conf, cu),1e-3);

 cu.reset();
 EXPECT_NEAR(173.149,of(conf, cu),1e-3);
 cu.change_rectangle(1,rects[2]);
 cu.remove_rectangle(0);
 EXPECT_NEAR(79.8879,of(conf, cu),1e-3);

 //  of.complete_operations();
//  EXPECT_NEAR(79.8879,of(),1e-3);
//  of.try_add_rectangle(rects[3]);
//  of.try_change_rectangle(0,rects[0]);
//  EXPECT_NEAR(167.785,of(),1e-3);
//  of.complete_operations();
//  EXPECT_NEAR(167.785,of(),1e-3);
}*/

MAKE_MAIN;
