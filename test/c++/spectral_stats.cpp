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

#include <utility>
#include <vector>

#include <triqs/mesh.hpp>

#include <nda/gtest_tools.hpp>

#include <som/spectral_stats.hpp>

using namespace nda;
using namespace som;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  if(mpi::has_env) {
    mpi::environment env(argc, argv);
    std::cout << "MPI environment detected" << std::endl;
    return RUN_ALL_TESTS();
  } else
    return RUN_ALL_TESTS();
}

namespace testing {
class spectral_stats_test : public ::testing::Test {
protected:
  som_core cont;
  triqs::mesh::refreq mesh = {-5.0, 5.0, 5};
  std::vector<std::pair<double, double>> interv = {{-5.0, -3.0},
                                                   {-1.0, 1.0},
                                                   {3.0, 5.0}};

  spectral_stats_test() {
    cont.data.resize(2, som_core::data_t{});

    cont.data[0].particular_solutions = {
        {configuration{{0, 1, 2}, {1, 4, 1}, {3, 2, 3}}, 0},
        {configuration{{-0.1, 1, 2}, {0.9, 4, 1}, {2.9, 2, 3}}, 0},
        {configuration{{-0.2, 1, 2}, {0.8, 4, 1}, {2.8, 2, 3}}, 0}};

    cont.data[1].particular_solutions = {
        {configuration{{0, 1, 2}, {-1, 4, 1}, {-3, 2, 3}}, 0},
        {configuration{{0.1, 1, 2}, {-0.9, 4, 1}, {-2.9, 2, 3}}, 0},
        {configuration{{0.2, 1, 2}, {-0.8, 4, 1}, {-2.8, 2, 3}}, 0}};
  }
};
} // namespace testing

using namespace testing;

TEST_F(spectral_stats_test, spectral_integral) {
  configuration c{{0, 1, 2}, {1, 4, 1}, {3, 2, 3}};

  EXPECT_NEAR(spectral_integral(2.0, 1.5, c, rectangle), 2.5, 1e-10);
  EXPECT_NEAR(spectral_integral(2.0, 1.5, c, lorentzian), 1.9842074284, 1e-10);
  EXPECT_NEAR(spectral_integral(2.0, 1.5, c, gaussian), 2.4419081106, 1e-10);

  triqs::arrays::vector<double> ref(mesh.size());

  ref = {0., 0., 1.7, 2.8, 0.3};
  EXPECT_ARRAY_NEAR(spectral_integral(mesh, c, rectangle), ref, 1e-10);
  ref = {0.1142532345, 0.3316005188, 1.3177758720, 1.8161452063, 0.6213209141};
  EXPECT_ARRAY_NEAR(spectral_integral(mesh, c, lorentzian), ref, 1e-10);
  ref = {0.0009945621, 0.2087447961, 1.5639704932, 2.3671263370, 0.6660794784};
  EXPECT_ARRAY_NEAR(spectral_integral(mesh, c, gaussian), ref, 1e-10);

  ref = {0., 2., 1.5};
  EXPECT_ARRAY_NEAR(spectral_integral(interv, c, rectangle), ref, 1e-10);
  ref = {0.1341663357, 1.4467315501, 1.2823782799};
  EXPECT_ARRAY_NEAR(spectral_integral(interv, c, lorentzian), ref, 1e-10);
  ref = {0.0018083638, 1.6740000753, 1.5908630342};
  EXPECT_ARRAY_NEAR(spectral_integral(interv, c, gaussian), ref, 1e-10);
}

TEST_F(spectral_stats_test, spectral_avg_mesh) {
  vector<double> ref(mesh.size());

  ref = {0., 0., 1.74, 2.88, 0.18};
  EXPECT_ARRAY_NEAR(spectral_avg(cont, 0, mesh, rectangle), ref, 1e-10);
  ref = {0.18, 2.88, 1.74, 0., 0.};
  EXPECT_ARRAY_NEAR(spectral_avg(cont, 1, mesh, rectangle), ref, 1e-10);

  ref = {0.1181354935, 0.3516508893, 1.3401479841, 1.8189434995, 0.5816390449};
  EXPECT_ARRAY_NEAR(spectral_avg(cont, 0, mesh, lorentzian), ref, 1e-10);
  ref = {0.5816390449, 1.8189434995, 1.3401479841, 0.3516508893, 0.1181354935};
  EXPECT_ARRAY_NEAR(spectral_avg(cont, 1, mesh, lorentzian), ref, 1e-10);

  ref = {0.0013523815, 0.2412393198, 1.6095326791, 2.3600231333, 0.5968717764};
  EXPECT_ARRAY_NEAR(spectral_avg(cont, 0, mesh, gaussian), ref, 1e-10);
  ref = {0.5968717764, 2.3600231333, 1.6095326791, 0.2412393198, 0.0013523815};
  EXPECT_ARRAY_NEAR(spectral_avg(cont, 1, mesh, gaussian), ref, 1e-10);
}

TEST_F(spectral_stats_test, spectral_avg_intervals) {
  vector<double> ref(interv.size());

  ref = {0., 2., 1.35};
  EXPECT_ARRAY_NEAR(spectral_avg(cont, 0, interv, rectangle), ref, 1e-10);
  ref = {1.35, 2., 0.};
  EXPECT_ARRAY_NEAR(spectral_avg(cont, 1, interv, rectangle), ref, 1e-10);

  ref = {0.1399390061, 1.4659822436, 1.1905339934};
  EXPECT_ARRAY_NEAR(spectral_avg(cont, 0, interv, lorentzian), ref, 1e-10);
  ref = {1.1905339934, 1.4659822436, 0.13993900613};
  EXPECT_ARRAY_NEAR(spectral_avg(cont, 1, interv, lorentzian), ref, 1e-10);

  ref = {0.0026140415807, 1.7088851639, 1.4631864139};
  EXPECT_ARRAY_NEAR(spectral_avg(cont, 0, interv, gaussian), ref, 1e-10);
  ref = {1.4631864139, 1.7088851639, 0.0026140415807};
  EXPECT_ARRAY_NEAR(spectral_avg(cont, 1, interv, gaussian), ref, 1e-10);
}

TEST_F(spectral_stats_test, spectral_disp_mesh) {
  vector<double> avg(mesh.size());
  vector<double> ref(mesh.size());

  ref = {0., 0., 0.0010666667, 0.0042666667, 0.0096};
  avg = spectral_avg(cont, 0, mesh, rectangle);
  EXPECT_ARRAY_NEAR(spectral_disp(cont, 0, mesh, avg, rectangle), ref, 1e-10);
  ref = {0.0096, 0.0042666667, 0.0010666667, 0., 0.};
  avg = spectral_avg(cont, 1, mesh, rectangle);
  EXPECT_ARRAY_NEAR(spectral_disp(cont, 1, mesh, avg, rectangle), ref, 1e-10);

  ref = {0.0000102231, 0.0002756819, 0.0003118767, 0.0000045915, 0.0010213268};
  avg = spectral_avg(cont, 0, mesh, lorentzian);
  EXPECT_ARRAY_NEAR(spectral_disp(cont, 0, mesh, avg, lorentzian), ref, 1e-10);
  ref = {0.0010213268, 0.0000045915, 0.0003118767, 0.0002756819, 0.0000102231};
  avg = spectral_avg(cont, 1, mesh, lorentzian);
  EXPECT_ARRAY_NEAR(spectral_disp(cont, 1, mesh, avg, lorentzian), ref, 1e-10);

  ref = {0.0000000937, 0.0007279478, 0.0013622646, 0.0000489348, 0.0031216267};
  avg = spectral_avg(cont, 0, mesh, gaussian);
  EXPECT_ARRAY_NEAR(spectral_disp(cont, 0, mesh, avg, gaussian), ref, 1e-10);
  ref = {0.0031216267, 0.0000489348, 0.0013622646, 0.0007279478, 0.0000000937};
  avg = spectral_avg(cont, 1, mesh, gaussian);
  EXPECT_ARRAY_NEAR(spectral_disp(cont, 1, mesh, avg, gaussian), ref, 1e-10);
}

TEST_F(spectral_stats_test, spectral_disp_intervals) {
  vector<double> avg(interv.size());
  vector<double> ref(interv.size());

  ref = {0., 0., 0.015};
  avg = spectral_avg(cont, 0, interv, rectangle);
  EXPECT_ARRAY_NEAR(spectral_disp(cont, 0, interv, avg, rectangle), ref, 1e-10);
  ref = {0.015, 0., 0.};
  avg = spectral_avg(cont, 1, interv, rectangle);
  EXPECT_ARRAY_NEAR(spectral_disp(cont, 1, interv, avg, rectangle), ref, 1e-10);

  ref = {0.0000227137, 0.0002172851, 0.0055832956};
  avg = spectral_avg(cont, 0, interv, lorentzian);
  EXPECT_ARRAY_NEAR(
      spectral_disp(cont, 0, interv, avg, lorentzian), ref, 1e-10);
  ref = {0.0055832956, 0.0002172851, 0.0000227138};
  avg = spectral_avg(cont, 1, interv, lorentzian);
  EXPECT_ARRAY_NEAR(
      spectral_disp(cont, 1, interv, avg, lorentzian), ref, 1e-10);

  ref = {0.0000004833, 0.0007690017, 0.0108434856};
  avg = spectral_avg(cont, 0, interv, gaussian);
  EXPECT_ARRAY_NEAR(spectral_disp(cont, 0, interv, avg, gaussian), ref, 1e-10);
  ref = {0.0108434856, 0.0007690017, 0.0000004833};
  avg = spectral_avg(cont, 1, interv, gaussian);
  EXPECT_ARRAY_NEAR(spectral_disp(cont, 1, interv, avg, gaussian), ref, 1e-10);
}

TEST_F(spectral_stats_test, spectral_corr_mesh) {
  vector<double> avg(mesh.size());
  array<double, 2> ref(mesh.size(), mesh.size());

  ref = {{0., 0., 0., 0., 0.},
         {0., 0., 0., 0., 0.},
         {0., 0., 0.0010666667, 0.0021333333, -0.0032},
         {0., 0., 0.0021333333, 0.0042666667, -0.0064},
         {0., 0., -0.0032, -0.0064, 0.0096}};
  avg = spectral_avg(cont, 0, mesh, rectangle);
  EXPECT_ARRAY_NEAR(spectral_corr(cont, 0, mesh, avg, rectangle), ref, 1e-10);
  ref = {{0.0096, -0.0064, -0.0032, 0., 0.},
         {-0.0064, 0.0042666667, 0.0021333333, 0., 0.},
         {-0.0032, 0.0021333333, 0.0010666667, 0., 0.},
         {0., 0., 0., 0., 0.},
         {0., 0., 0., 0., 0.}};
  avg = spectral_avg(cont, 1, mesh, rectangle);
  EXPECT_ARRAY_NEAR(spectral_corr(cont, 1, mesh, avg, rectangle), ref, 1e-10);

  ref = {
      {0.0000102231, 0.0000530856, 0.0000562950, 0.0000040839, -0.0001021032},
      {0.0000530856, 0.0002756819, 0.0002921184, 0.0000209473, -0.0005300060},
      {0.0000562950, 0.0002921185, 0.0003118767, 0.0000248466, -0.0005639645},
      {0.0000040839, 0.0000209473, 0.0000248466, 0.0000045915, -0.0000429414},
      {-0.0001021032,
       -0.0005300060,
       -0.00056396445954,
       -0.0000429414,
       0.0010213268}};
  avg = spectral_avg(cont, 0, mesh, lorentzian);
  EXPECT_ARRAY_NEAR(spectral_corr(cont, 0, mesh, avg, lorentzian), ref, 1e-10);
  ref = {
      {0.0010213268,
       -0.0000429414,
       -0.0005639645,
       -0.0005300060,
       -0.0001021032},
      {-0.0000429414, 0.0000045915, 0.0000248466, 0.00002094731, 0.0000040839},
      {-0.0005639645, 0.0000248466, 0.0003118767, 0.00029211849, 0.0000562950},
      {-0.0005300060, 0.0000209473, 0.0002921185, 0.00027568189, 0.0000530856},
      {-0.0001021032, 0.0000040839, 0.0000562950, 0.00005308557, 0.0000102231}};
  avg = spectral_avg(cont, 1, mesh, lorentzian);
  EXPECT_ARRAY_NEAR(spectral_corr(cont, 1, mesh, avg, lorentzian), ref, 1e-10);

  ref = {
      {0.0000000937, 0.0000082505, 0.0000112547, -0.0000021098, -0.0000170273},
      {0.0000082505, 0.0007279478, 0.0009949408, -0.0001842657, -0.0015056906},
      {0.0000112547, 0.0009949408, 0.0013622646, -0.0002495023, -0.0020621152},
      {-0.0000021098, -0.0001842657, -0.0002495023, 0.0000489348, 0.0003770601},
      {-0.0000170273,
       -0.0015056906,
       -0.0020621152,
       0.0003770601,
       0.0031216267}};
  avg = spectral_avg(cont, 0, mesh, gaussian);
  EXPECT_ARRAY_NEAR(spectral_corr(cont, 0, mesh, avg, gaussian), ref, 1e-10);
  ref = {
      {0.0031216267, 0.0003770601, -0.0020621152, -0.0015056906, -0.0000170273},
      {0.0003770601, 0.0000489348, -0.0002495023, -0.0001842657, -0.0000021098},
      {-0.0020621152, -0.0002495023, 0.0013622646, 0.0009949408, 0.0000112547},
      {-0.0015056906, -0.0001842657, 0.0009949408, 0.0007279478, 0.0000082505},
      {-0.0000170273, -0.0000021098, 0.0000112547, 0.0000082505, 0.0000000937}};
  avg = spectral_avg(cont, 1, mesh, gaussian);
  EXPECT_ARRAY_NEAR(spectral_corr(cont, 1, mesh, avg, gaussian), ref, 1e-10);
}

TEST_F(spectral_stats_test, spectral_corr_intervals) {
  vector<double> avg(interv.size());
  array<double, 2> ref(interv.size(), interv.size());

  ref = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.015}};
  avg = spectral_avg(cont, 0, interv, rectangle);
  EXPECT_ARRAY_NEAR(spectral_corr(cont, 0, interv, avg, rectangle), ref, 1e-10);
  ref = {{0.015, 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
  avg = spectral_avg(cont, 1, interv, rectangle);
  EXPECT_ARRAY_NEAR(spectral_corr(cont, 1, interv, avg, rectangle), ref, 1e-10);

  ref = {{0.0000227138, 0.0000694781, -0.0003560028},
         {0.0000694781, 0.0002172851, -0.0010930463},
         {-0.0003560028, -0.0010930463, 0.0055832956}};
  avg = spectral_avg(cont, 0, interv, lorentzian);
  EXPECT_ARRAY_NEAR(
      spectral_corr(cont, 0, interv, avg, lorentzian), ref, 1e-10);
  ref = {{0.0055832956, -0.0010930463, -0.0003560028},
         {-0.0010930463, 0.0002172851, 0.0000694781},
         {-0.0003560028, 0.0000694781, 0.0000227138}};
  avg = spectral_avg(cont, 1, interv, lorentzian);
  EXPECT_ARRAY_NEAR(
      spectral_corr(cont, 1, interv, avg, lorentzian), ref, 1e-10);

  ref = {{0.0000004833, 0.0000191015, -0.0000721109},
         {0.0000191015, 0.0007690017, -0.0028844581},
         {-0.0000721109, -0.0028844581, 0.0108434856}};
  avg = spectral_avg(cont, 0, interv, gaussian);
  EXPECT_ARRAY_NEAR(spectral_corr(cont, 0, interv, avg, gaussian), ref, 1e-10);
  ref = {{0.0108434856, -0.0028844581, -0.0000721109},
         {-0.0028844581, 0.0007690017, 0.0000191015},
         {-0.0000721109, 0.0000191015, 0.0000004833}};
  avg = spectral_avg(cont, 1, interv, gaussian);
  EXPECT_ARRAY_NEAR(spectral_corr(cont, 1, interv, avg, gaussian), ref, 1e-10);
}
