#pragma once

#include <vector>
#include <cmath>
#include <utility>

#include <triqs/arrays/vector.hpp>

#include "configuration.hpp"
#include "config_update.hpp"

namespace som {

template<typename KernelType>
class objective_function {

 using rhs_type = typename KernelType::result_type;
 using mesh_type = typename KernelType::mesh_type;

 // Integral kernel
 KernelType kern;
 // The right-hand side of the Fredholm integral equation
 rhs_type const& rhs;
 // Error bars of the RHS
 rhs_type const& error_bars;

public:

 objective_function(mesh_type const& mesh,
                    rhs_type const& rhs,
                    rhs_type const& error_bars) :
  kern(mesh),
  rhs(rhs), error_bars(error_bars) {}

 double operator()(configuration const& c) const {
  return sum(abs((kern(c) - rhs) / error_bars));
 }

 double operator()(configuration const& c, config_update const& cu) const {
  return sum(abs((kern(c, cu) - rhs) / error_bars));
 }

 KernelType const& get_kernel() { return kern; }
};

}
