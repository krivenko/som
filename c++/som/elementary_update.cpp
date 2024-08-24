/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#ifdef EXT_DEBUG
#include <cmath>
#include <iostream>

#include <triqs/utility/numeric_ops.hpp>
#endif

#include "elementary_update.hpp"
#include "solution_worker.hpp"

namespace som {

using namespace triqs::mc_tools;

template <typename KernelType>
std::pair<bool, double>
elementary_update<KernelType>::optimize_parameter_change(double dxi,
                                                         double dxi_min,
                                                         double dxi_max) {
  new_objf_value[full] = data.objf(update[full]);
  new_objf_value[half] = data.objf(update[half]);

  double a =
      2 *
      (data.temp_objf_value - 2 * new_objf_value[half] + new_objf_value[full]) /
      (dxi * dxi);
  if(a <= 0)
    return std::make_pair(false, 0); // dxi_opt is going to be a maximum
  double b = (4 * new_objf_value[half] - 3 * data.temp_objf_value -
              new_objf_value[full]) /
             dxi;

  double dxi_opt = -b / (2 * a);

  // dxi_opt is outside the definition domain
  if(dxi_opt < dxi_min || dxi_opt >= dxi_max)
    return std::make_pair(false, 0);

  else
    return std::make_pair(true, dxi_opt);
}

// Make the final choice between :math:`\delta\xi`, :math:`\delta\xi/2`,
// and, optionally, :math:`\delta\xi_{opt}`.
template <typename KernelType>
void elementary_update<KernelType>::select_parameter_change(bool consider_opt) {
  if(consider_opt) {
    new_objf_value[opt] = data.objf(update[opt]);
    selected_parameter_change = opt;
    if(new_objf_value[full] < new_objf_value[opt])
      selected_parameter_change = full;
    if(new_objf_value[half] < new_objf_value[full])
      selected_parameter_change = half;
  } else {
    selected_parameter_change =
        (new_objf_value[full] < new_objf_value[half]) ? full : half;
  }
}

// Manually set selected_parameter_change and update new_objf_value
template <typename KernelType>
void elementary_update<KernelType>::set_parameter_change(
    parameter_change_t parameter_change) {
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
  new_objf_value[parameter_change] = data.objf(update[parameter_change]);
  selected_parameter_change = parameter_change;
}

// Returns Metropolis ratio for the currently selected update
template <typename KernelType>
double elementary_update<KernelType>::transition_probability() const {
  double old_d = data.temp_objf_value;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
  double new_d = new_objf_value[selected_parameter_change];
  return new_d < old_d ? 1.0 : data.Z(old_d / new_d);
}

template <typename KernelType>
#ifdef EXT_DEBUG
elementary_update<KernelType>::elementary_update(
    mc_data<KernelType>& data,
    random_generator& rng,
    cache_index& ci,
    std::shared_ptr<cc_update_t> cc_update,
    std::pair<double, double> energy_window,
    double width_min,
    double weight_min)
   : data(data)
   , rng(rng)
   , ci(ci)
   , kern(data.objf.get_kernel())
   , energy_window(std::move(energy_window))
   , width_min(width_min)
   , weight_min(weight_min)
#else
elementary_update<KernelType>::elementary_update(
    mc_data<KernelType>& data,
    random_generator& rng,
    cache_index& ci,
    std::shared_ptr<cc_update_t> cc_update)
   : data(data)
   , rng(rng)
   , ci(ci)
   , kern(data.objf.get_kernel())
#endif
   , dist(rng, data.gamma)
   , cc_update(std::move(cc_update))
   , update{config_update(data.temp_conf, ci),
            config_update(data.temp_conf, ci),
            config_update(data.temp_conf, ci)}
   , new_objf_value{0, 0, 0} {
}

//----------------

template <typename KernelType> double elementary_update<KernelType>::accept() {
  // Update temporary configuration and its objective function value
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
  update[selected_parameter_change].apply();
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
  kern.cache_copy(update[selected_parameter_change], data.temp_conf);
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
  data.temp_objf_value = new_objf_value[selected_parameter_change];

  // Reset all considered updates
  update[full].reset();
  update[half].reset();
  update[opt].reset();

#ifdef EXT_DEBUG
  std::cerr << "* Elementary update accepted\n";
  std::cerr << "Temporary configuration: size = " << data.temp_conf.size()
            << ", norm = " << data.temp_conf.norm()
            << ", χ = " << std::sqrt(data.temp_objf_value) << '\n';
  using triqs::utility::is_zero;
  if(!is_zero(data.temp_conf.norm() - data.global_conf.norm(), 1e-10))
    TRIQS_RUNTIME_ERROR
        << "Elementary update has changed norm of the solution!";

  for(auto& r : data.temp_conf) {
    if(r.norm() < weight_min)
      TRIQS_RUNTIME_ERROR << "Rectangle is too small: " << r;
    if(r.width < width_min)
      TRIQS_RUNTIME_ERROR << "Rectangle is too narrow: " << r;
    if(r.height < 0)
      TRIQS_RUNTIME_ERROR << "Rectangle with negative height: " << r;
    if(r.left() < energy_window.first || r.right() > energy_window.second)
      TRIQS_RUNTIME_ERROR << "Rectangle is not within energy window: " << r;
  }
#endif

  // Copy this temporary configuration to the globally selected configuration
  // if its objective function value is smaller
  if(data.temp_objf_value < data.global_objf_value) {
#ifdef EXT_DEBUG
    std::cerr << "Copying temporary configuration to global configuration "
              << "(χ(temp) = " << std::sqrt(data.temp_objf_value)
              << ", χ(global) = " << std::sqrt(data.global_objf_value) << ")\n";
#endif
    data.global_conf = data.temp_conf;
    kern.cache_copy(data.temp_conf, data.global_conf);
    data.global_objf_value = data.temp_objf_value;
  }

#ifdef EXT_DEBUG
  std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
#endif

  // Update internal state of the distribution function
  ++data.Z;

  return 1;
}

//----------------
template <typename KernelType> void elementary_update<KernelType>::reject() {

#ifdef EXT_DEBUG
  std::cerr << "* Elementary update rejected\n";
  std::cerr << "Temporary configuration: size = " << data.temp_conf.size()
            << ", norm = " << data.temp_conf.norm()
            << ", χ = " << std::sqrt(data.temp_objf_value) << '\n';
  std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
#endif

  // Reset all considered updates
  update[full].reset();
  update[half].reset();
  update[opt].reset();

  // Update internal state of the distribution function
  ++data.Z;
}

INSTANTIATE_CLASS_FOR_EACH_KERNEL(elementary_update)

} // namespace som
