#pragma once
#include <triqs/parameters.hpp>
#include <triqs/gfs.hpp>
#include <triqs/arrays.hpp>

using namespace triqs::gfs;
using namespace triqs::utility;
using namespace triqs::arrays;

class my_class {

 double Beta;
 block_gf<imtime> g1;
 gf<refreq> g2;
 gf<imtime> g3;
 array<double,1> Arr;
 array<double,1> ar;
 matrix<double> m1;

 public:
 my_class(triqs::utility::parameters p);

 // specify all required and optional parameters and generate help from them
 triqs::utility::parameter_defaults constructor_defaults() const;

 // Access to some data, for Python ...
 // do NOT add const here : python use a non const object and non const view
 
 double beta() { return Beta; }
 block_gf_view<imtime> g1_view() { return g1; }
 gf_view<refreq> g2_view() { return g2; }
 gf_view<imtime> g3_view() { return g3; }
 array_view<double,1> Arr_view() { return Arr; }
 array_view<double,1> ar_view() { return ar; }
 matrix_view<double> m1_view() { return m1; }
};

