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
#include <triqs/test_tools/arrays.hpp>

#include <triqs/gfs.hpp>

#include <som/configuration.hpp>

using namespace triqs::arrays;
using namespace triqs::gfs;
using namespace som;

cache_index ci;

configuration
    conf1({{-2.0, 2.6, 0.3, ci}, {1.3, 2.6, 0.6, ci}, {2.0, 2.6, 0.7, ci}}, ci);
std::string conf1_str =
    "(c:-2, w:2.6, h:0.3),(c:1.3, w:2.6, h:0.6),(c:2, w:2.6, h:0.7)";
configuration
    conf2({{-3.0, 2.6, 0.3, ci}, {-1.6, 2.4, 0.1, ci}, {2.8, 2.8, 0.55, ci}},
          ci);
std::string conf2_str =
    "(c:-3, w:2.6, h:0.3),(c:-1.6, w:2.4, h:0.1),(c:2.8, w:2.8, h:0.55)";

TEST(configuration, Evaluation) {
  gf_mesh<refreq> mesh(-5.0, 5.0, 11);
  vector<double> ref1 = {.0, .0, 0.3, 0.3, 0.3, 0.6, 1.3, 1.3, 0.7, .0, .0};
  vector<double> ref2 = {.0, 0.3, 0.3, 0.4, 0.1, .0, .0, 0.55, 0.55, 0.55, .0};

  for(int k = 0; k < mesh.size(); ++k) {
    EXPECT_CLOSE(conf1(mesh[k]), ref1(k));
    EXPECT_CLOSE(conf2(mesh[k]), ref2(k));
  }
}

TEST(configuration, Print) {
  ASSERT_PRINT(conf1_str, conf1);
  ASSERT_PRINT(conf2_str, conf2);
}

TEST(configuration, Detached) {
  configuration conf_det({{-2.0, 2.6, 0.3}, {1.3, 2.6, 0.6}, {2.0, 2.6, 0.7}});
  ASSERT_PRINT(conf1_str, conf_det);
}

TEST(configuration, Arithmetics) {
  EXPECT_PRINT(conf1_str + "," + conf2_str, conf1 + conf2);

  EXPECT_PRINT("(c:-2, w:2.6, h:0.6),(c:1.3, w:2.6, h:1.2),(c:2, w:2.6, h:1.4)",
               2.0 * conf1);
  EXPECT_PRINT(
      "(c:-3, w:2.6, h:0.6),(c:-1.6, w:2.4, h:0.2),(c:2.8, w:2.8, h:1.1)",
      conf2 * 2.0);

  conf1 += conf2;
  EXPECT_PRINT(conf1_str + "," + conf2_str, conf1);
  conf2 *= 3.0;
  EXPECT_PRINT(
      "(c:-3, w:2.6, h:0.9),(c:-1.6, w:2.4, h:0.3),(c:2.8, w:2.8, h:1.65)",
      conf2);
}

TEST(configuration, normalize) {
  conf2.normalize(2.56);
  ASSERT_PRINT(conf2_str, conf2);
}

MAKE_MAIN
