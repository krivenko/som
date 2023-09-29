/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2023 Igor Krivenko <igor.s.krivenko@gmail.com>
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
#include <cmath>
#include <vector>

// clang-format off
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>
// clang-format on

#include <triqs/mesh.hpp>

#include <som/kernels/fermiongf_imfreq.hpp>
#include <som/kernels/fermiongf_imtime.hpp>
#include <som/solution_functionals/objective_function.hpp>

using namespace nda;
using namespace som;

class objective_function_test : public ::testing::Test {

protected:
  const double beta = 2;

  cache_index ci;
  std::vector<rectangle> rects = {{-2, 2.6, 0.3, ci},
                                  {1.3, 2.6, 0.6, ci},
                                  {-0.5, 2.6, 0.5, ci},
                                  {2, 2.6, 0.7, ci}};

public:
  objective_function_test() = default;
};

#define MAKE_TEST_CHANGE(TEST_CASE, TOL)                                       \
  TEST_F(TEST_CASE, Change) {                                                  \
    configuration conf({rects[0], rects[2]}, ci);                              \
                                                                               \
    EXPECT_NEAR(chi2_ref_R13, of(conf), TOL);                                  \
                                                                               \
    config_update cu(conf, ci);                                                \
    cu.change_rectangle(1, rects[3]);                                          \
    EXPECT_NEAR(chi2_ref_R14, of(cu), TOL);                                    \
                                                                               \
    cu.reset();                                                                \
    EXPECT_NEAR(chi2_ref_R13, of(cu), TOL);                                    \
    cu.change_rectangle(0, rects[1]);                                          \
    EXPECT_NEAR(chi2_ref_R23, of(cu), TOL);                                    \
                                                                               \
    cu.apply();                                                                \
    of.get_kernel().cache_copy(cu, conf);                                      \
    EXPECT_NEAR(chi2_ref_R23, of(conf), TOL);                                  \
  }

#define MAKE_TEST_ADD(TEST_CASE, TOL)                                          \
  TEST_F(TEST_CASE, Add) {                                                     \
    configuration conf({rects[0], rects[1]}, ci);                              \
                                                                               \
    EXPECT_NEAR(chi2_ref_R12, of(conf), TOL);                                  \
                                                                               \
    config_update cu(conf, ci);                                                \
    cu.add_rectangle(rects[2]);                                                \
    EXPECT_NEAR(chi2_ref_R123, of(cu), TOL);                                   \
                                                                               \
    cu.reset();                                                                \
    EXPECT_NEAR(chi2_ref_R12, of(cu), TOL);                                    \
    cu.add_rectangle(rects[3]);                                                \
    EXPECT_NEAR(chi2_ref_R124, of(cu), TOL);                                   \
                                                                               \
    cu.apply();                                                                \
    of.get_kernel().cache_copy(cu, conf);                                      \
    EXPECT_NEAR(chi2_ref_R124, of(conf), TOL);                                 \
  }

#define MAKE_TEST_REMOVE(TEST_CASE, TOL)                                       \
  TEST_F(TEST_CASE, Remove) {                                                  \
    configuration conf({rects[0], rects[1], rects[2]}, ci);                    \
                                                                               \
    EXPECT_NEAR(chi2_ref_R123, of(conf), TOL);                                 \
                                                                               \
    config_update cu(conf, ci);                                                \
    cu.remove_rectangle(1);                                                    \
    EXPECT_NEAR(chi2_ref_R13, of(cu), TOL);                                    \
                                                                               \
    cu.reset();                                                                \
    EXPECT_NEAR(chi2_ref_R123, of(cu), TOL);                                   \
    cu.remove_rectangle(2);                                                    \
    EXPECT_NEAR(chi2_ref_R12, of(cu), TOL);                                    \
                                                                               \
    cu.apply();                                                                \
    of.get_kernel().cache_copy(cu, conf);                                      \
    EXPECT_NEAR(chi2_ref_R12, of(conf), TOL);                                  \
  }

#define MAKE_TEST_MULTIPLE(TEST_CASE, TOL)                                     \
  TEST_F(TEST_CASE, Multiple) {                                                \
    configuration conf({rects[0], rects[1]}, ci);                              \
                                                                               \
    EXPECT_NEAR(chi2_ref_R12, of(conf), TOL);                                  \
                                                                               \
    config_update cu(conf, ci);                                                \
    cu.add_rectangle(rects[3]);                                                \
    cu.change_rectangle(1, rects[2]);                                          \
    EXPECT_NEAR(chi2_ref_R134, of(cu), TOL);                                   \
                                                                               \
    cu.reset();                                                                \
    EXPECT_NEAR(chi2_ref_R12, of(cu), TOL);                                    \
    cu.change_rectangle(1, rects[2]);                                          \
    cu.remove_rectangle(0);                                                    \
    EXPECT_NEAR(chi2_ref_R3, of(cu), TOL);                                     \
                                                                               \
    cu.apply();                                                                \
    of.get_kernel().cache_copy(cu, conf);                                      \
    EXPECT_NEAR(chi2_ref_R3, of(conf), TOL);                                   \
    cu.add_rectangle(rects[3]);                                                \
    cu.change_rectangle(0, rects[0]);                                          \
    EXPECT_NEAR(chi2_ref_R14, of(cu), TOL);                                    \
    cu.apply();                                                                \
    EXPECT_NEAR(chi2_ref_R14, of(conf), TOL);                                  \
  }

//
// Imaginary time
//

class objective_function_imtime_test : public objective_function_test {
protected:
  triqs::mesh::imtime mesh;
  kernel<FermionGf, triqs::mesh::imtime> kern;
  array<double, 1> g;

public:
  objective_function_imtime_test()
     : mesh(beta, triqs::mesh::Fermion, 11), kern(mesh), g(mesh.size()) {
    for(auto pt : mesh) {
      size_t i = pt.data_index();
      auto tau = double(pt);
      g(i) = -0.5 * (std::exp(-tau * 1.3) / (1 + std::exp(-beta * 1.3)) +
                     std::exp(tau * 0.7) / (1 + std::exp(beta * 0.7)));
    }
  }
};

//
// Error bars
//

class objective_function_imtime_error_bars_test
   : public objective_function_imtime_test {
  array<double, 1> make_error_bars() {
    array<double, 1> res(mesh.size());
    for(auto pt : mesh) {
      size_t i = pt.data_index();
      auto tau = double(pt);
      res(i) = 0.05 * std::exp(-std::abs(tau - 0.5 * beta));
    }
    return res;
  }

  array<double, 1> s;

protected:
  objective_function<kernel<FermionGf, triqs::mesh::imtime>> of;

  double const chi2_ref_R3 = 94.5020111;
  double const chi2_ref_R12 = 393.205112;
  double const chi2_ref_R13 = 586.079405;
  double const chi2_ref_R23 = 991.805377;
  double const chi2_ref_R14 = 577.132629;
  double const chi2_ref_R123 = 1715.98578;
  double const chi2_ref_R124 = 2797.05285;
  double const chi2_ref_R134 = 1860.57671;
  double const chi2_ref_R24 = 2532.14361;

public:
  objective_function_imtime_error_bars_test()
     : s(make_error_bars()), of(kern, g, s) {}
};

MAKE_TEST_CHANGE(objective_function_imtime_error_bars_test, 1e-5)
MAKE_TEST_ADD(objective_function_imtime_error_bars_test, 1e-5)
MAKE_TEST_REMOVE(objective_function_imtime_error_bars_test, 1e-5)
MAKE_TEST_MULTIPLE(objective_function_imtime_error_bars_test, 1e-5)

//
// Covariance matrix
//

class objective_function_imtime_cov_matrix_test
   : public objective_function_imtime_test {
  matrix<double> make_cov_matrix() {
    auto res = matrix<double>::zeros({mesh.size(), mesh.size()});
    for(auto i = 0; i < mesh.size(); ++i) res(i, i) = 0.1;
    for(auto i = 0; i < mesh.size() - 1; ++i) {
      res(i, i + 1) = 0.2;
      res(i + 1, i) = 0.2;
    }
    return res;
  }

  matrix<double> cov_matrix;
  double l = 0.1;

protected:
  objective_function<kernel<FermionGf, triqs::mesh::imtime>> of;

  double const chi2_ref_R3 = 0.230858697;
  double const chi2_ref_R12 = 0.797454077;
  double const chi2_ref_R13 = 1.31186285;
  double const chi2_ref_R23 = 2.29158832;
  double const chi2_ref_R14 = 1.12347611;
  double const chi2_ref_R123 = 3.80930694;
  double const chi2_ref_R124 = 5.80561801;
  double const chi2_ref_R134 = 3.87333242;
  double const chi2_ref_R24 = 5.43485068;

public:
  objective_function_imtime_cov_matrix_test()
     : cov_matrix(make_cov_matrix()), of(kern, g, cov_matrix, l) {}
};

MAKE_TEST_CHANGE(objective_function_imtime_cov_matrix_test, 1e-7)
MAKE_TEST_ADD(objective_function_imtime_cov_matrix_test, 1e-7)
MAKE_TEST_REMOVE(objective_function_imtime_cov_matrix_test, 1e-7)
MAKE_TEST_MULTIPLE(objective_function_imtime_cov_matrix_test, 1e-7)

//
// Imaginary frequency
//

class objective_function_imfreq_test : public objective_function_test {
protected:
  triqs::mesh::imfreq mesh;
  kernel<FermionGf, triqs::mesh::imfreq> kern;
  array<dcomplex, 1> g;

public:
  objective_function_imfreq_test()
     : mesh(beta,
            triqs::mesh::Fermion,
            11,
            triqs::mesh::imfreq::option::positive_frequencies_only)
     , kern(mesh)
     , g(mesh.size()) {
    for(auto pt : mesh) {
      size_t i = pt.data_index();
      auto w = dcomplex(pt);
      g(i) = 0.5 * (1.0 / (w - 1.3) + 1.0 / (w + 0.7));
    }
  }
};

//
// Error bars
//

class objective_function_imfreq_error_bars_test
   : public objective_function_imfreq_test {
  array<dcomplex, 1> make_error_bars() {
    array<dcomplex, 1> res(mesh.size());
    for(auto pt : mesh) {
      size_t i = pt.data_index();
      auto w = dcomplex(pt);
      res(i) = 0.05 / std::abs(w);
    }
    return res;
  }

  array<dcomplex, 1> s;

protected:
  objective_function<kernel<FermionGf, triqs::mesh::imfreq>> of;

  double const chi2_ref_R3 = 42.3978067;
  double const chi2_ref_R12 = 626.145646;
  double const chi2_ref_R13 = 451.906742;
  double const chi2_ref_R23 = 1290.22582;
  double const chi2_ref_R14 = 865.954239;
  double const chi2_ref_R123 = 2532.73609;
  double const chi2_ref_R124 = 3544.25872;
  double const chi2_ref_R134 = 2962.41627;
  double const chi2_ref_R24 = 2109.81232;

public:
  objective_function_imfreq_error_bars_test()
     : s(make_error_bars()), of(kern, g, s) {}
};

MAKE_TEST_CHANGE(objective_function_imfreq_error_bars_test, 1e-5)
MAKE_TEST_ADD(objective_function_imfreq_error_bars_test, 1e-5)
MAKE_TEST_REMOVE(objective_function_imfreq_error_bars_test, 1e-5)
MAKE_TEST_MULTIPLE(objective_function_imfreq_error_bars_test, 1e-5)

//
// Covariance matrix
//

class objective_function_imfreq_cov_matrix_test
   : public objective_function_imfreq_test {
  matrix<dcomplex> make_cov_matrix() {
    auto res = matrix<dcomplex>::zeros({mesh.size(), mesh.size()});
    using namespace std::complex_literals;
    for(auto i = 0; i < mesh.size(); ++i) res(i, i) = 0.1;
    for(auto i = 0; i < mesh.size() - 1; ++i) {
      res(i, i + 1) = -0.2i;
      res(i + 1, i) = 0.2i;
    }
    return res;
  }

  matrix<dcomplex> cov_matrix;
  double l = 0.1;

protected:
  objective_function<kernel<FermionGf, triqs::mesh::imfreq>> of;

  double const chi2_ref_R3 = 0.0510053912;
  double const chi2_ref_R12 = 0.236197806;
  double const chi2_ref_R13 = 0.261843449;
  double const chi2_ref_R23 = 0.753306531;
  double const chi2_ref_R14 = 0.263218328;
  double const chi2_ref_R123 = 1.19940677;
  double const chi2_ref_R124 = 1.58761383;
  double const chi2_ref_R134 = 1.14429089;
  double const chi2_ref_R24 = 1.30382031;

public:
  objective_function_imfreq_cov_matrix_test()
     : cov_matrix(make_cov_matrix()), of(kern, g, cov_matrix, l) {}
};

MAKE_TEST_CHANGE(objective_function_imfreq_cov_matrix_test, 1e-8)
MAKE_TEST_ADD(objective_function_imfreq_cov_matrix_test, 1e-8)
MAKE_TEST_REMOVE(objective_function_imfreq_cov_matrix_test, 1e-8)
MAKE_TEST_MULTIPLE(objective_function_imfreq_cov_matrix_test, 1e-8)
