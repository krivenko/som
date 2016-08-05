/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016 by I. Krivenko
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
#include <triqs/arrays.hpp>

#include "kernels/fermiongf_imtime.hpp"
#include "kernels/fermiongf_imfreq.hpp"
#include "kernels/fermiongf_legendre.hpp"

using namespace som;
using namespace triqs::gfs;
using triqs::arrays::vector;

cache_index ci;

TEST(FermionGf, imtime) {

 h5::file file("fermiongf_imtime.h5", H5F_ACC_RDWR);

 configuration conf(ci);
 h5_read(file, "rects", conf, ci);
 vector<double> beta;
 h5_read(file, "beta", beta);

 for(double b : beta) {
  array<double,2> results;
  h5_read(file, "results/beta" + std::to_string(int(b)), results);

  gf_mesh<imtime> mesh(b, Fermion, second_dim(results));
  kernel<FermionGf,imtime> kern(mesh);
  ci.invalidate_all();

  for(int i : range(conf.size()))
   EXPECT_TRUE(array_are_close(results(i,range()), kern(conf[i]), 1e-10));
 }
}

TEST(FermionGf, imfreq) {

 h5::file file("fermiongf_imfreq.h5", H5F_ACC_RDWR);

 configuration conf(ci);
 h5_read(file, "rects", conf, ci);
 vector<double> beta;
 h5_read(file, "beta", beta);

 for(double b : beta) {
  array<double,2> results_re, results_im;
  h5_read(file, "results/beta" + std::to_string(int(b)) + "/re", results_re);
  h5_read(file, "results/beta" + std::to_string(int(b)) + "/im", results_im);
  array<dcomplex,2> results = results_re + 1_j*results_im;

  gf_mesh<imfreq> mesh(b, Fermion, second_dim(results));
  kernel<FermionGf,imfreq> kern(mesh);
  ci.invalidate_all();

  for(int i : range(conf.size()))
   EXPECT_TRUE(array_are_close(results(i,range()), kern(conf[i]), 1e-10));
 }
}

TEST(FermionGf, legendre) {

 h5::file file("fermiongf_legendre.h5", H5F_ACC_RDWR);

 configuration conf(ci);
 h5_read(file, "rects", conf, ci);
 vector<double> beta;
 h5_read(file, "beta", beta);

 for(double b : beta) {
  array<double,2> results;
  h5_read(file, "results/beta" + std::to_string(int(b)), results);

  gf_mesh<legendre> mesh(b, Fermion, second_dim(results));
  kernel<FermionGf,legendre> kern(mesh);
  ci.invalidate_all();

  for(int i : range(conf.size()))
   EXPECT_TRUE(array_are_close(results(i,range()), kern(conf[i]), 1e-10));
 }
}

MAKE_MAIN;
