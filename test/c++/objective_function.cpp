/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
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
#include <vector>

#include <triqs/test_tools/arrays.hpp>
#include <triqs/gfs.hpp>

#include "kernels/fermiongf_imtime.hpp"
#include "objective_function.hpp"

using namespace som;
using namespace triqs::arrays;
using namespace triqs::gfs;

const double beta = 2;
cache_index ci;
std::vector<rectangle> rects {{-2,2.6,0.3,ci},
                              {1.3,2.6,0.6,ci},
                              {-0.5,2.6,0.5,ci},
                              {2,2.6,0.7,ci}};

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

auto g = GF();
auto s = S();

kernel<FermionGf,imtime> kern(mesh);

objective_function<kernel<FermionGf,imtime>> of(kern, g, s);

TEST(objective_function,Change) {
 configuration conf({rects[0],rects[2]},ci);

 EXPECT_NEAR(172.784,of(conf),1e-3);

 config_update cu(conf,ci);
 cu.change_rectangle(1, rects[3]);
 EXPECT_NEAR(167.785,of(cu),1e-3);

 cu.reset();
 EXPECT_NEAR(172.784,of(cu),1e-3);
 cu.change_rectangle(0,rects[1]);
 EXPECT_NEAR(301.098,of(cu),1e-3);

 cu.apply();
 of.get_kernel().cache_copy(cu, conf);
 EXPECT_NEAR(301.098,of(conf),1e-3);
}

TEST(objective_function,Add) {
 configuration conf({rects[0],rects[1]},ci);

 EXPECT_NEAR(173.149,of(conf),1e-3);

 config_update cu(conf,ci);
 cu.add_rectangle(rects[2]);
 EXPECT_NEAR(400.832,of(cu),1e-3);

 cu.reset();
 EXPECT_NEAR(173.149,of(cu),1e-3);
 cu.add_rectangle(rects[3]);
 EXPECT_NEAR(405.862,of(cu),1e-3);

 cu.apply();
 of.get_kernel().cache_copy(cu, conf);
 EXPECT_NEAR(405.862,of(conf),1e-3);
}

TEST(objective_function,Remove) {
 configuration conf({rects[0],rects[1],rects[2]},ci);

 EXPECT_NEAR(400.832,of(conf),1e-3);

 config_update cu(conf,ci);
 cu.remove_rectangle(1);
 EXPECT_NEAR(172.784,of(cu),1e-3);

 cu.reset();
 EXPECT_NEAR(400.832,of(cu),1e-3);
 cu.remove_rectangle(2);
 EXPECT_NEAR(173.149,of(cu),1e-3);

 cu.apply();
 of.get_kernel().cache_copy(cu, conf);
 EXPECT_NEAR(173.149,of(conf),1e-3);
}


TEST(objective_function,Multiple) {
 configuration conf({rects[0],rects[1]},ci);

 EXPECT_NEAR(173.149,of(conf),1e-3);

 config_update cu(conf,ci);
 cu.add_rectangle(rects[3]);
 cu.change_rectangle(1,rects[2]);
 EXPECT_NEAR(395.468,of(cu),1e-3);

 cu.reset();
 EXPECT_NEAR(173.149,of(cu),1e-3);
 cu.change_rectangle(1,rects[2]);
 cu.remove_rectangle(0);
 EXPECT_NEAR(79.8879,of(cu),1e-3);

 cu.apply();
 of.get_kernel().cache_copy(cu, conf);
 EXPECT_NEAR(79.8879,of(conf),1e-3);
 cu.add_rectangle(rects[3]);
 cu.change_rectangle(0,rects[0]);
 EXPECT_NEAR(167.785,of(cu),1e-3);
 cu.apply();
 EXPECT_NEAR(167.785,of(conf),1e-3);
}

MAKE_MAIN;
