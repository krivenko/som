/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2026 Igor Krivenko
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

#include <complex>
#include <string>

// clang-format off
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>
// clang-format on

#include <h5/h5.hpp>

#include <triqs/utility/is_complex.hpp>

using namespace std::complex_literals;
using namespace triqs::mesh;
using nda::array;
using nda::range;
using nda::vector;
using triqs::is_complex;

using namespace som;

namespace som {

// h5_read_mathematica_array(), complex arrays
template <typename ArrayType>
ArrayType h5_read_mathematica_array(h5::group g, std::string const& name)
  requires(is_complex<typename ArrayType::value_type>::value)
{
  array<typename ArrayType::value_type, ArrayType::rank> res_re, res_im;
  h5_read(g, name + "/re", res_re);
  h5_read(g, name + "/im", res_im);
  return res_re + 1i * res_im;
}

// h5_read_mathematica_array(), real arrays
template <typename ArrayType>
ArrayType h5_read_mathematica_array(h5::group g, std::string const& name)
  requires(!is_complex<typename ArrayType::value_type>::value)
{
  ArrayType res;
  h5_read(g, name, res);
  return res;
}

// test_kernel()
template <typename KernelType>
void test_kernel(std::string const& filename,
                 cache_index& ci,
                 double tolerance) {

  h5::file file(filename, 'r');

  configuration conf(ci);
  h5_read(file, "rects", conf);
  vector<double> beta;
  h5_read(file, "beta", beta);

  for(double b : beta) {
    using mesh_t = typename KernelType::mesh_type;
    using scalar_t = typename KernelType::result_type::value_type;

    std::string group_name = "results/beta" + std::to_string(int(b));
    auto results =
        h5_read_mathematica_array<array<scalar_t, 2>>(file, group_name);
    mesh_t mesh(b,
                KernelType::kind == ZeroTemp
                    ? triqs::gfs::Fermion
                    : observable_statistics(KernelType::kind),
                second_dim(results));
    KernelType kern(mesh);
    ci.invalidate_all();

    for(auto i : range(long(conf.size())))
      EXPECT_TRUE(
          array_are_close(results(i, range::all), kern(conf[i]), tolerance));
  }
}

} // namespace som
