#include <vector>

#include <test_tools.hpp>
#include <triqs/gfs.hpp>

#include "kernels.hpp"
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
  res(i) = -0.5*(exp(-tau*1.3)/(1 + exp(-beta*1.3)) +
                 exp(tau*0.7)/(1 + exp(beta*0.7)));
 }
 return res;
}

array<double,1> S() {
 array<double,1> res(mesh.size());
 for(auto pt : mesh) {
  int i = pt.index();
  double tau = double(pt);
  res(i) = 0.05*exp(-abs(tau - 0.5*beta));
 }
 return res;
}

using obj_function = objective_function<kernel<FermionicGf,imtime>>;

TEST(objective_function,Change) {
 configuration conf {rects[0],rects[2]};
 auto gf = GF();
 auto s = S();

 obj_function of(mesh, gf, s, conf);

 std::cout << of() << std::endl;
}

TEST(objective_function,Add) {
}

TEST(objective_function,Remove) {
}

TEST(objective_function,Multiple) {
}

MAKE_MAIN;
