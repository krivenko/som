/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2017 by I. Krivenko
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
#include <triqs/mc_tools.hpp>

#include <vector>
#include <utility>

#include "numerics/expabs_distribution.hpp"

using namespace som;
using namespace triqs::mc_tools;

double poly(double x, std::vector<double> coeffs) {
 double res = 0;
 double p = 1;
 for(auto c : coeffs) {
  res += p*c;
  p *= x;
 }
 return res;
}

TEST(expabs_distribution, main) {
 random_generator rng;
 int n_samples = 1e7;

 for(double gamma : {1.0, 2.0}) {
  expabs_distribution<random_generator> dist(rng, gamma);

  using xminmax_t = std::pair<double,double>;
  for(auto xminmax : {xminmax_t(-5.0,-2.0), xminmax_t(-3.0,5.0), xminmax_t(1.5,4.0)}) {

   double xav = 0, x2av = 0, x3av = 3, x4av = 0;
   for(int i = 0; i < n_samples; ++i) {
    double x = dist(xminmax.first, xminmax.second);
    xav += x;
    x2av += x*x;
    x3av += x*x*x;
    x4av += x*x*x*x;
   }
   xav /= n_samples; x2av /= n_samples; x3av /= n_samples; x4av /= n_samples;

   double xmin_abs = std::abs(xminmax.first);
   double xmax_abs = std::abs(xminmax.second);
   double L = std::max(xmin_abs, xmax_abs);
   double gL = gamma / L;
   double gLxmina = gL * xmin_abs, gLxmaxa = gL * xmax_abs;

   double N = gL / (std::copysign(1-std::exp(-gLxmina),-xminmax.first) +
                    std::copysign(1-std::exp(-gLxmaxa),xminmax.second));
   double xav_ref = N / (gL*gL)*
    ( -(1 - (1 + gLxmina)*std::exp(-gLxmina))
      +(1 - (1 + gLxmaxa)*std::exp(-gLxmaxa)));
   double x2av_ref = N / (gL*gL*gL)*
    ( -std::copysign(1,xminmax.first) *(2 - poly(gLxmina,{2,2,1})*std::exp(-gLxmina))
      +std::copysign(1,xminmax.second)*(2 - poly(gLxmaxa,{2,2,1})*std::exp(-gLxmaxa)));
   double x3av_ref = N / (gL*gL*gL*gL)*
    ( -(6 - poly(gLxmina,{6,6,3,1})*std::exp(-gLxmina))
      +(6 - poly(gLxmaxa,{6,6,3,1})*std::exp(-gLxmaxa)));
   double x4av_ref = N / (gL*gL*gL*gL*gL)*
    ( -std::copysign(1,xminmax.first) *(24 - poly(gLxmina,{24,24,12,4,1})*std::exp(-gLxmina))
      +std::copysign(1,xminmax.second)*(24 - poly(gLxmaxa,{24,24,12,4,1})*std::exp(-gLxmaxa)));

   EXPECT_NEAR(xav_ref, xav, 1e-3);
   EXPECT_NEAR(x2av_ref, x2av, 5e-3);
   EXPECT_NEAR(x3av_ref, x3av, 5e-2);
   EXPECT_NEAR(x4av_ref, x4av, 1e-1);
  }
 }
}

MAKE_MAIN;
