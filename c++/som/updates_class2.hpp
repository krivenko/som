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

// Adding a new rectangle (elementary update D)
template <typename KernelType>
class update_insert : public elementary_update<KernelType> {

  std::pair<double, double> energy_window;
  double width_min;
  double weight_min;
  double cnew_min, cnew_max;
  int max_rects;

  using eu = elementary_update<KernelType>;

public:
  update_insert(mc_data<KernelType>& data,
                triqs::mc_tools::random_generator& rng, cache_index& ci,
                std::pair<double, double> energy_window, double width_min,
                double weight_min, int max_rects)
     :
#ifdef EXT_DEBUG
     elementary_update<KernelType>(data, rng, ci, energy_window, width_min,
                                   weight_min)
     ,
#else
     elementary_update<KernelType>(data, rng, ci)
     ,
#endif
     energy_window(std::move(energy_window))
     , width_min(width_min)
     , weight_min(weight_min)
     , cnew_min(energy_window.first + width_min / 2)
     , cnew_max(energy_window.second - width_min / 2)
     , max_rects(max_rects) {
  }

  double attempt();
};

//------------------------------------------------------------------------------

// Removing one rectangle and shifting another (elementary update E)
template <typename KernelType>
class update_remove_shift : public elementary_update<KernelType> {

  std::pair<double, double> energy_window;

  using eu = elementary_update<KernelType>;

public:
  update_remove_shift(mc_data<KernelType>& data,
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

// Splitting a rectangle and shifting the pieces (elementary update F)
template <typename KernelType>
class update_split_shift : public elementary_update<KernelType> {

  std::pair<double, double> energy_window;
  double width_min;
  double weight_min;
  int max_rects;

  using eu = elementary_update<KernelType>;

public:
  update_split_shift(mc_data<KernelType>& data,
                     triqs::mc_tools::random_generator& rng, cache_index& ci,
                     std::pair<double, double> energy_window, double width_min,
                     double weight_min, int max_rects)
     :
#ifdef EXT_DEBUG
     elementary_update<KernelType>(data, rng, ci, energy_window, width_min,
                                   weight_min)
     ,
#else
     elementary_update<KernelType>(data, rng, ci)
     ,
#endif
     energy_window(std::move(energy_window))
     , width_min(width_min)
     , weight_min(weight_min)
     , max_rects(max_rects) {
  }

  double attempt();
};

//------------------------------------------------------------------------------

// Glueing rectangles (elementary update G)
template <typename KernelType>
class update_glue_shift : public elementary_update<KernelType> {

  std::pair<double, double> energy_window;

  using eu = elementary_update<KernelType>;

public:
  update_glue_shift(mc_data<KernelType>& data,
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

EXTERN_TEMPLATE_CLASS_FOR_EACH_KERNEL(update_insert);
EXTERN_TEMPLATE_CLASS_FOR_EACH_KERNEL(update_remove_shift);
EXTERN_TEMPLATE_CLASS_FOR_EACH_KERNEL(update_split_shift);
EXTERN_TEMPLATE_CLASS_FOR_EACH_KERNEL(update_glue_shift);

} // namespace som
