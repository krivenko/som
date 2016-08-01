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

#include "kernels/bosoncorrsym_imfreq.hpp"

using namespace som;
using namespace triqs::gfs;
using triqs::arrays::vector;

cache_index ci;
configuration conf({{1.30000001, 2.6, 0.3,ci},
                    {2.7, 2.6, 0.4,ci},
                    {3.5, 2.6, 0.5,ci},
                    {5.3, 2.6, 0.6,ci},
                    {6.0, 2.6, 0.7,ci}},ci);

TEST(BosonCorrSym, imfreq) {

 double beta = 1;
 gf_mesh<imfreq> mesh(beta,Boson,5);
 kernel<BosonCorrSym,imfreq> kern(mesh);
 ci.invalidate_all();

 std::vector<vector<dcomplex>> ref = {
  {0.4965634224467141,0.025745261602540856,0.006909053409302071,0.003113715066794663,0.0017601271165260985},
  {0.6620845632622854,0.10580309418218788,0.03098652585148805,0.014235964352031234,0.008103960398589181},
  {0.8276057040778567,0.1964752722001335,0.06138600338250142,0.028643275964017088,0.016400184502754774},
  {0.9931268448934282,0.4097515593832911,0.1509977093598627,0.07371801378520859,0.042957411554051944},
  {1.158647985708999,0.5485863422424838,0.2158177138361782,0.10755230455854847,0.06319925550919453}
 };

 for(int i = 0; i < conf.size(); ++i) EXPECT_TRUE(array_are_close(ref[i],kern(conf[i]),1e-6));
}

MAKE_MAIN;
