#include <vector>
#include <triqs/test_tools/arrays.hpp>

#include "kernels.hpp"

using namespace som;
using namespace triqs::gfs;
using triqs::arrays::vector;

TEST(Kernels,FermionicGf_imtime) {

 double beta = 1;
 configuration conf = {{-2.0,2.6,0.3},
                       {-1.3,2.6,0.4},
                       {-0.5,2.6,0.5},
                       {1.3,2.6,0.6},
                       {2.0,2.6,0.7}};

 std::vector<vector<double>> ref = {
 {-0.11009,-0.12894,-0.15171,-0.17935,-0.21309,-0.25448,-0.30549,-0.36872,-0.44746,-0.54600,-0.66991},
 {-0.24860,-0.27281,-0.30080,-0.33328,-0.37115,-0.41546,-0.46754,-0.52898,-0.60174,-0.68822,-0.79140},
 {-0.50906,-0.51974,-0.53328,-0.54995,-0.57010,-0.59411,-0.62244,-0.65562,-0.69427,-0.73910,-0.79094},
 {-1.18710,-1.03234,-0.90261,-0.79347,-0.70131,-0.62319,-0.55672,-0.49992,-0.45119,-0.40922,-0.37290},
 {-1.56312,-1.27400,-1.04408,-0.86036,-0.71283,-0.59378,-0.49721,-0.41849,-0.35399,-0.30086,-0.25688}};

 gf_mesh<imtime> mesh(beta,Fermion,11);

 kernel<FermionicGf,imtime> kern(mesh);

 for(int i = 0; i < conf.size(); ++i)
  EXPECT_TRUE(array_are_close(ref[i],kern(conf[i]),1e-5));

}

MAKE_MAIN;
