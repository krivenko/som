/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
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
#include <sstream>

// clang-format off
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>
// clang-format on

#include <triqs/mesh.hpp>

#include <som/configuration.hpp>

using namespace nda;
using namespace som;

#define EXPECT_PRINT(X, Y)                                                     \
  {                                                                            \
    std::stringstream ss;                                                      \
    ss << Y;                                                                   \
    EXPECT_EQ(X, ss.str());                                                    \
  }

class configuration_test : public ::testing::Test {
protected:
  cache_index ci;
  configuration conf1 = {
      {{-2.0, 2.6, 0.3, ci}, {1.3, 2.6, 0.6, ci}, {2.0, 2.6, 0.7, ci}},
      ci};
  configuration conf2 = {
      {{-3.0, 2.6, 0.3, ci}, {-1.6, 2.4, 0.1, ci}, {2.8, 2.8, 0.55, ci}},
      ci};

  std::string const conf1_str =
      "(c:-2, w:2.6, h:0.3),(c:1.3, w:2.6, h:0.6),(c:2, w:2.6, h:0.7)";
  std::string const conf2_str =
      "(c:-3, w:2.6, h:0.3),(c:-1.6, w:2.4, h:0.1),(c:2.8, w:2.8, h:0.55)";

public:
  configuration_test() = default;
};

TEST_F(configuration_test, Evaluation) {
  triqs::mesh::refreq mesh(-5.0, 5.0, 11);
  vector<double> ref1 = {.0, .0, 0.3, 0.3, 0.3, 0.6, 1.3, 1.3, 0.7, .0, .0};
  vector<double> ref2 = {.0, 0.3, 0.3, 0.4, 0.1, .0, .0, 0.55, 0.55, 0.55, .0};

  for(int k = 0; k < mesh.size(); ++k) {
    EXPECT_CLOSE(conf1(mesh[k]), ref1(k));
    EXPECT_CLOSE(conf2(mesh[k]), ref2(k));
  }
}

TEST_F(configuration_test, Print) {
  EXPECT_PRINT(conf1_str, conf1);
  EXPECT_PRINT(conf2_str, conf2);
}

TEST_F(configuration_test, Detached) {
  configuration conf_det({{-2.0, 2.6, 0.3}, {1.3, 2.6, 0.6}, {2.0, 2.6, 0.7}});
  EXPECT_PRINT(conf1_str, conf_det);
}

TEST_F(configuration_test, Arithmetics) {
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

TEST_F(configuration_test, normalize) {
  conf2.normalize(2.56);
  EXPECT_PRINT(conf2_str, conf2);
}

TEST_F(configuration_test, prune) {
  configuration conf({{-2.0, 2.6, 0.3, ci},
                      {-1.0, 1e-5, 1e5, ci},
                      {1.3, 2.6, 0.6, ci},
                      {2.0, 2.6, 0.7, ci},
                      {3.0, 0.01, 0.01, ci}},
                     ci);
  conf.prune(1e-4, 1e-3);
  configuration ref(
      {{-2.0, 2.6, 0.3, ci}, {1.3, 2.6, 0.6, ci}, {2.0, 2.6, 0.7, ci}});
  EXPECT_EQ(conf, ref);
}

TEST_F(configuration_test, redistribute_small_rects_weight) {
  configuration conf1({{1.0, 0.1, 0, ci},
                       {2.0, 0.2, 0.01, ci},
                       {3.0, 0.3, 2.0, ci},
                       {4.0, 0.4, -0.01, ci},
                       {5.0, 0.5, 0.015, ci}},
                      ci);
  conf1.redistribute_small_rects_weight(0.08);
  EXPECT_EQ(conf1, configuration({{3.0, 0.3, 2.0183333333333335}}));

  configuration conf2({{1.0, 0.1, 0, ci},
                       {2.0, 0.2, 0, ci},
                       {3.0, 0.3, 2, ci},
                       {4.0, 0.4, 0, ci},
                       {5.0, 0.5, 0, ci}},
                      ci);
  conf2.redistribute_small_rects_weight(0.001);
  EXPECT_EQ(conf2, configuration({{3.0, 0.3, 2.0}}));

  configuration conf3({{1.0, 0.1, 0, ci},
                       {2.0, 0.2, 0, ci},
                       {3.0, 0.3, 0, ci},
                       {4.0, 0.4, 0, ci},
                       {5.0, 0.5, 0, ci}},
                      ci);
  conf3.redistribute_small_rects_weight(0.001);
  EXPECT_EQ(conf3, configuration());
}

TEST_F(configuration_test, make_nonoverlapping) {
  configuration conf({{0.5, 6.0, 1.0}, {3.1, 1.8, 2.0}, {-0.5, 1.0, 3.0}}, ci);
  configuration ref({{-2.75, 0.5, 0.0, ci},
                     {-1.75, 1.5, 1.0, ci},
                     {-0.5, 1.0, 4.0, ci},
                     {1.1, 2.2, 1.0, ci},
                     {2.85, 1.3, 3.0, ci},
                     {3.75, 0.5, 2.0, ci},
                     {4.5, 1.0, 0.0, ci}},
                    ci);
  auto nonoverlapping = make_nonoverlapping(conf, {-3.0, 5.0});
  EXPECT_EQ(nonoverlapping, ref);

  configuration ref_width_min({{-1.75, 1.5, 1.0, ci},
                               {-0.5, 1.0, 4.0, ci},
                               {1.1, 2.2, 1.0, ci},
                               {2.85, 1.3, 3.0, ci},
                               {4.5, 1.0, 0.0, ci}},
                              ci);
  auto nonoverlapping_width_min = make_nonoverlapping(conf, {-3.0, 5.0}, 0.6);
  EXPECT_EQ(nonoverlapping_width_min, ref_width_min);
}

TEST_F(configuration_test, strip_rect_heights) {
  configuration conf(
      {{0.5, 6.0, 1.0, ci}, {3.1, 1.8, 2.0, ci}, {-0.5, 1.0, 3.0, ci}}, ci);

  vector<double> heights(3);
  conf.strip_rect_heights(heights);

  configuration ref(
      {{0.5, 6.0, 1.0, ci}, {3.1, 1.8, 1.0, ci}, {-0.5, 1.0, 1.0, ci}}, ci);
  vector<double> heights_ref = {1.0, 2.0, 3.0};
  EXPECT_EQ(conf, ref);
  EXPECT_EQ(heights, heights_ref);
}

TEST_F(configuration_test, update_rect_heights) {
  configuration conf(
      {{0.5, 6.0, 1.0, ci}, {3.1, 1.8, 1.0, ci}, {-0.5, 1.0, 1.0, ci}}, ci);

  vector<double> heights = {1.0, 2.0, 3.0};
  conf.update_rect_heights(heights);

  configuration ref(
      {{0.5, 6.0, 1.0, ci}, {3.1, 1.8, 2.0, ci}, {-0.5, 1.0, 3.0, ci}}, ci);
  EXPECT_EQ(conf, ref);
}
