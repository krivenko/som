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

#include <nda/gtest_tools.hpp>

#include <triqs/mesh.hpp>

#include <som/kernels/fermiongf_imtime.hpp>
#include <som/solution_functionals/objective_function.hpp>

using namespace nda;
using namespace som;

class objective_function_test : public ::testing::Test {

  const double beta = 2;
  triqs::mesh::imtime mesh;
  kernel<FermionGf, triqs::mesh::imtime> kern;

  array<double, 1> GF() {
    array<double, 1> res(mesh.size());
    for(auto pt : mesh) {
      size_t i = pt.linear_index();
      auto tau = double(pt);
      res(i) = -0.5 * (std::exp(-tau * 1.3) / (1 + std::exp(-beta * 1.3)) +
                       std::exp(tau * 0.7) / (1 + std::exp(beta * 0.7)));
    }
    return res;
  }

  array<double, 1> S() {
    array<double, 1> res(mesh.size());
    for(auto pt : mesh) {
      size_t i = pt.linear_index();
      auto tau = double(pt);
      res(i) = 0.05 * std::exp(-std::abs(tau - 0.5 * beta));
    }
    return res;
  }

  array<double, 1> g;
  array<double, 1> s;

protected:
  cache_index ci;
  std::vector<rectangle> rects = {{-2, 2.6, 0.3, ci},
                                  {1.3, 2.6, 0.6, ci},
                                  {-0.5, 2.6, 0.5, ci},
                                  {2, 2.6, 0.7, ci}};

  objective_function<kernel<FermionGf, triqs::mesh::imtime>> of;

public:
  objective_function_test()
     : mesh(beta, triqs::mesh::Fermion, 11)
     , kern(mesh)
     , g(GF())
     , s(S())
     , of(kern, g, s) {}
};

TEST_F(objective_function_test, Change) {
  configuration conf({rects[0], rects[2]}, ci);

  EXPECT_NEAR(586.079, of(conf), 1e-3);

  config_update cu(conf, ci);
  cu.change_rectangle(1, rects[3]);
  EXPECT_NEAR(577.133, of(cu), 1e-3);

  cu.reset();
  EXPECT_NEAR(586.079, of(cu), 1e-3);
  cu.change_rectangle(0, rects[1]);
  EXPECT_NEAR(991.805, of(cu), 1e-3);

  cu.apply();
  of.get_kernel().cache_copy(cu, conf);
  EXPECT_NEAR(991.805, of(conf), 1e-3);
}

TEST_F(objective_function_test, Add) {
  configuration conf({rects[0], rects[1]}, ci);

  EXPECT_NEAR(393.205, of(conf), 1e-3);

  config_update cu(conf, ci);
  cu.add_rectangle(rects[2]);
  EXPECT_NEAR(1715.986, of(cu), 1e-3);

  cu.reset();
  EXPECT_NEAR(393.205, of(cu), 1e-3);
  cu.add_rectangle(rects[3]);
  EXPECT_NEAR(2797.053, of(cu), 1e-3);

  cu.apply();
  of.get_kernel().cache_copy(cu, conf);
  EXPECT_NEAR(2797.053, of(conf), 1e-3);
}

TEST_F(objective_function_test, Remove) {
  configuration conf({rects[0], rects[1], rects[2]}, ci);

  EXPECT_NEAR(1715.986, of(conf), 1e-3);

  config_update cu(conf, ci);
  cu.remove_rectangle(1);
  EXPECT_NEAR(586.079, of(cu), 1e-3);

  cu.reset();
  EXPECT_NEAR(1715.986, of(cu), 1e-3);
  cu.remove_rectangle(2);
  EXPECT_NEAR(393.205, of(cu), 1e-3);

  cu.apply();
  of.get_kernel().cache_copy(cu, conf);
  EXPECT_NEAR(393.205, of(conf), 1e-3);
}

TEST_F(objective_function_test, Multiple) {
  configuration conf({rects[0], rects[1]}, ci);

  EXPECT_NEAR(393.205, of(conf), 1e-3);

  config_update cu(conf, ci);
  cu.add_rectangle(rects[3]);
  cu.change_rectangle(1, rects[2]);
  EXPECT_NEAR(1860.577, of(cu), 1e-3);

  cu.reset();
  EXPECT_NEAR(393.205, of(cu), 1e-3);
  cu.change_rectangle(1, rects[2]);
  cu.remove_rectangle(0);
  EXPECT_NEAR(94.502, of(cu), 1e-3);

  cu.apply();
  of.get_kernel().cache_copy(cu, conf);
  EXPECT_NEAR(94.502, of(conf), 1e-3);
  cu.add_rectangle(rects[3]);
  cu.change_rectangle(0, rects[0]);
  EXPECT_NEAR(577.133, of(cu), 1e-3);
  cu.apply();
  EXPECT_NEAR(577.133, of(conf), 1e-3);
}
