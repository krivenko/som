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
#include <vector>
#include <triqs/test_tools/arrays.hpp>

#include "kernels/bosoncorr_imfreq.hpp"

using namespace som;
using namespace triqs::gfs;
using triqs::arrays::vector;

cache_index ci;
configuration conf({{1.30000001, 2.6, 0.3,ci},
                    {2.7, 2.6, 0.4,ci},
                    {3.5, 2.6, 0.5,ci},
                    {5.3, 2.6, 0.6,ci},
                    {6.0, 2.6, 0.7,ci}},ci);

TEST(BosonCorr, imfreq) {

 double beta = 1;
 gf_mesh<imfreq> mesh(beta,Boson,5);
 kernel<BosonCorr,imfreq> kern(mesh);
 ci.invalidate_all();

 std::vector<vector<dcomplex>> ref = {
  {0.248282+0_j,0.0128726+0.0474171_j,0.00345453+0.0251504_j,0.00155686+0.0169624_j,0.000880064+0.0127742_j},
  {0.331042+0_j,0.0529015+0.116714_j,0.0154933+0.0673401_j,0.00711798+0.0462551_j,0.00405198+0.0350666_j},
  {0.413803+0_j,0.0982376+0.17203_j,0.030693+0.106003_j,0.0143216+0.0739495_j,0.00820009+0.0563859_j},
  {0.496563+0_j,0.204876+0.241986_j,0.0754989+0.17652_j,0.036859+0.128878_j,0.0214787+0.100015_j},
  {0.579324+0_j,0.274293+0.28697_j,0.107909+0.223784_j,0.0537762+0.166807_j,0.0315996+0.130539_j}
 };

 for(int i = 0; i < conf.size(); ++i) EXPECT_TRUE(array_are_close(ref[i],kern(conf[i]),1e-6));
}

MAKE_MAIN;
