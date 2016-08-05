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
                    {-2.0,2.6,0.3,ci},
                    {-1.3,2.6,0.4,ci},
                    {-0.5,2.6,0.5,ci},
                    {1.3,2.6,0.6,ci},
                    {2.0,2.6,0.7,ci}},ci);

TEST(BosonCorr, imfreq) {

 double beta = 1;
 gf_mesh<imfreq> mesh(beta,Boson,5);
 kernel<BosonCorr,imfreq> kern(mesh);
 ci.invalidate_all();

 std::vector<vector<dcomplex>> ref = {
  {0.248282+0_j,0.0128726+0.0474171_j,0.00345453+0.0251504_j,0.00155686+0.0169624_j,0.000880064+0.0127742_j},
  {0.248282+0_j,0.0246847-0.0693824_j,0.00688917-0.0381538_j,0.00313087-0.0259301_j,0.00177519-0.0195817_j},
  {0.331042+0_j,0.0171635-0.0632227_j,0.00460604-0.0335338_j,0.00207581-0.0226166_j,0.00117342-0.0170323_j},
  {0.413803+0_j,0.00815389-0.0313973_j,0.0021071-0.0162652_j,0.000942423-0.0109169_j,0.000531294-0.00820715_j},
  {0.496563+0_j,0.0257453+0.0948341_j,0.00690905+0.0503007_j,0.00311372+0.0339248_j,0.00176013+0.0255485_j},
  {0.579324+0_j,0.0575976+0.161892_j,0.0160747+0.0890255_j,0.00730535+0.0605035_j,0.00414211+0.0456906_j}
 };

 for(int i = 0; i < conf.size(); ++i) EXPECT_TRUE(array_are_close(ref[i],kern(conf[i]),1e-6));
}

MAKE_MAIN;