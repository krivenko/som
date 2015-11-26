#pragma once

#include <vector>
#include <cmath>
#include <utility>

#include <triqs/arrays/vector.hpp>

#include "configuration.hpp"

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

 // Deviation vector, see eq. (30),
 // and the value it will have after a call to complete_operation()
 rhs_type deviation;

 // Value of the function and the value it will have
 double value;

 // Partial LHS: the kernel applied to each rectangle
 std::vector<rhs_type> partial_lhs;

 // List of changed partial_lhs with their indices
 // 0 <= index < partial_lhs.size(): replace partial_lhs at index
 // index == partial_lhs.size(): add partial_lhs
 // index < 0: remove partial_lhs at -index.
 std::vector<std::pair<int,rhs_type>> changed_partial_lhs;

 // Values of 'deviation' and 'value' after a call to complete_operation()
 mutable rhs_type new_deviation;
 mutable double new_value;

public:

 objective_function(mesh_type const& mesh,
                    rhs_type const& rhs,
                    rhs_type const& error_bars,
                    configuration const& init_config = {}
                    ) :
  kern(mesh),
  rhs(rhs), error_bars(error_bars),
  deviation(mesh.size()),
  new_deviation(mesh.size()) {

  deviation() = new_deviation() = 0;
  changed_partial_lhs.reserve(2);

  partial_lhs.reserve(init_config.max_size());
  for(auto const& r : init_config) {
   partial_lhs.push_back(kern(r));
   deviation += partial_lhs.back();
  }
  deviation -= rhs;
  deviation /= error_bars;

  value = sum(abs(deviation));
 }

 void try_change_rectangle(int i, rectangle const& r) {
  changed_partial_lhs.emplace_back(i,kern(r));
 }
 void try_add_rectangle(rectangle const& r) {
  changed_partial_lhs.emplace_back(partial_lhs.size(),kern(r));
 }
 void try_remove_rectangle(int i) {
  changed_partial_lhs.emplace_back(-i-1,{});
 }

 double operator()() const {
  if(changed_partial_lhs.size()) {
   new_deviation() = 0;
   for(auto const& plhs : changed_partial_lhs) {
    int index = plhs.first;
    if(index == partial_lhs.size()) { // add rectangle
     new_deviation += plhs.second;
    } else if(index < 0) { // remove rectangle
     new_deviation -= partial_lhs[-index-1];
    } else { // change rectangle
     new_deviation += plhs.second - partial_lhs[index];
    }
   }
   new_deviation = new_deviation / error_bars + deviation;
   new_value = sum(abs(new_deviation));
   return new_value;
  } else
   return value;
 }

 void cancel_operations() {
  changed_partial_lhs.resize(0);
 }
 void complete_operation() {
  for(auto const& plhs : changed_partial_lhs) {
   int index = plhs.first;
   if(index == partial_lhs.size()) {
    partial_lhs.push_back(plhs.second);
   } else if(index < 0) {
    std::swap(partial_lhs[-index-1],partial_lhs.back());
    partial_lhs.pop_back();
   } else {
    partial_lhs[index] = plhs.second;
   }
  }
  changed_partial_lhs.resize(0);
  value = new_value;
  deviation() = new_deviation();
 }

};

}
