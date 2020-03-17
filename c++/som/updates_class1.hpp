/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#include "elementary_update.hpp"

namespace som {

// Shift of rectangle (elementary update A)
template <typename KernelType>
class update_shift : public elementary_update<KernelType> {

  std::pair<double, double> energy_window;

  using eu = elementary_update<KernelType>;

public:
  update_shift(mc_data<KernelType>& data,
               triqs::mc_tools::random_generator& rng, cache_index& ci,
               std::pair<double, double> energy_window,
#ifdef EXT_DEBUG
               double width_min, double weight_min)
     : elementary_update<KernelType>(data, rng, ci, energy_window, width_min,
                                     weight_min)
     ,
#else
               double, double)
     : elementary_update<KernelType>(data, rng, ci)
     ,
#endif
     energy_window(std::move(energy_window)) {
  }

  double attempt();
};

//------------------------------------------------------------------------------

// Change of width without change of weight (elementary update B)
template <typename KernelType>
class update_change_width : public elementary_update<KernelType> {

  std::pair<double, double> energy_window;
  double width_min;

  using eu = elementary_update<KernelType>;

public:
  update_change_width(mc_data<KernelType>& data,
                      triqs::mc_tools::random_generator& rng, cache_index& ci,
                      std::pair<double, double> energy_window,
#ifdef EXT_DEBUG
                      double width_min, double weight_min)
     : elementary_update<KernelType>(data, rng, ci, energy_window, width_min,
                                     weight_min)
     ,
#else
                      double width_min, double)
     : elementary_update<KernelType>(data, rng, ci)
     ,
#endif
     energy_window(std::move(energy_window))
     , width_min(width_min) {
  }

  double attempt();
};

//------------------------------------------------------------------------------

// Change of weight of two rectangles (elementary update C)
template <typename KernelType>
class update_change_weight2 : public elementary_update<KernelType> {

  double weight_min;

  using eu = elementary_update<KernelType>;

public:
  update_change_weight2(mc_data<KernelType>& data,
                        triqs::mc_tools::random_generator& rng, cache_index& ci,
#ifdef EXT_DEBUG
                        std::pair<double, double> const& energy_window,
                        double width_min, double weight_min)
     : elementary_update<KernelType>(data, rng, ci, std::move(energy_window),
                                     width_min, weight_min)
     ,
#else
                        std::pair<double, double> const&, double,
                        double weight_min)
     : elementary_update<KernelType>(data, rng, ci)
     ,
#endif
     weight_min(weight_min) {
  }

  double attempt();
};

EXTERN_TEMPLATE_CLASS_FOR_EACH_KERNEL(update_shift);
EXTERN_TEMPLATE_CLASS_FOR_EACH_KERNEL(update_change_width);
EXTERN_TEMPLATE_CLASS_FOR_EACH_KERNEL(update_change_weight2);

} // namespace som
