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
#pragma once

#include <string>
#include <triqs/h5.hpp>
#include <triqs/test_tools/arrays.hpp>
#include <triqs/utility/is_complex.hpp>
#include <triqs/utility/c14.hpp>

using namespace som;
using triqs::arrays::vector;
using triqs::arrays::array;
using triqs::is_complex;

namespace som {

// h5_read_mathematica_array(), complex arrays
template<typename ArrayType>
std14::enable_if_t<is_complex<typename ArrayType::value_type>::value, ArrayType>
h5_read_mathematica_array(h5::group g, std::string const& name) {
 array<typename ArrayType::value_type, ArrayType::rank> res_re, res_im;
 h5_read(g, name + "/re", res_re);
 h5_read(g, name + "/im", res_im);
 return res_re + 1_j*res_im;
}

// h5_read_mathematica_array(), real arrays
template<typename ArrayType>
std14::enable_if_t<!is_complex<typename ArrayType::value_type>::value, ArrayType>
h5_read_mathematica_array(h5::group g, std::string const& name) {
 ArrayType res;
 h5_read(g, name, res);
 return res;
}

// test_kernel()
template<typename KernelType>
void test_kernel(std::string const& filename, cache_index & ci, double tolerance) {

 h5::file file(filename, H5F_ACC_RDONLY);

 configuration conf(ci);
 h5_read(file, "rects", conf, ci);
 vector<double> beta;
 h5_read(file, "beta", beta);

 for(double b : beta) {
  using mesh_t = typename KernelType::mesh_type;
  using scalar_t = typename KernelType::result_type::value_type;

  std::string group_name = "results/beta" + std::to_string(int(b));
  auto results = h5_read_mathematica_array<array<scalar_t,2>>(file, group_name);
  mesh_t mesh(b,
              KernelType::kind == ZeroTemp ? Fermion : observable_statistics(KernelType::kind),
              second_dim(results));
  KernelType kern(mesh);
  ci.invalidate_all();

  for(int i : range(conf.size()))
   EXPECT_TRUE(array_are_close(results(i,range()), kern(conf[i]), tolerance));
 }
}

}
