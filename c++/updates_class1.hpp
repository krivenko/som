/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
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
template<typename KernelType> class update_shift : public elementary_update<KernelType> {

 std::pair<double,double> energy_window;

 using eu = elementary_update<KernelType>;

public:

 update_shift(mc_data<KernelType> & data, random_generator & rng, cache_index & ci,
              std::pair<double,double> energy_window, double width_min, double weight_min) :
  INIT_EU_BASE(data,rng,ci,energy_window,width_min,weight_min),
  energy_window(energy_window)
 {}

 double attempt() {

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "* Proposing update_shift (A)" << std::endl;
#endif

  int t = eu::rng(eu::data.temp_conf.size());
  auto const& rect = eu::data.temp_conf[t];

  double dc_min = energy_window.first + rect.width/2 - rect.center;
  double dc_max = energy_window.second - rect.width/2 - rect.center;
  double dc = eu::generate_parameter_change(dc_min, dc_max);
  if(dc == 0) return 0;

  eu::update[eu::full].change_rectangle(t, {rect.center + dc, rect.width, rect.height, eu::ci});
  eu::update[eu::half].change_rectangle(t, {rect.center + dc/2, rect.width, rect.height, eu::ci});

  auto dc_opt = eu::optimize_parameter_change(dc, dc_min, dc_max);

#ifdef EXT_DEBUG
  std::cerr << "Selected rectangle: " << rect << " [" << t << "]" << std::endl;
  std::cerr << "dc_min = " << dc_min << ", dc_max = " << dc_max
            << ", dc = " << dc << ", dc_opt = " << dc_opt.second << std::endl;
#endif

  if(dc_opt.first) {
   dc = dc_opt.second;
   eu::update[eu::opt].change_rectangle(t, {rect.center + dc, rect.width, rect.height, eu::ci});
  }

  eu::select_parameter_change(dc_opt.first);

#ifdef EXT_DEBUG
  std::cerr << "selected_parameter_change = " << eu::selected_parameter_change << std::endl;
#endif

  return eu::transition_probability();
 }

};

//-------------------------------------------------------------------------------------

// Change of width without change of weight (elementary update B)
template<typename KernelType> class update_change_width : public elementary_update<KernelType> {

 std::pair<double,double> energy_window;
 double width_min;

 using eu = elementary_update<KernelType>;

public:

 update_change_width(mc_data<KernelType> & data, random_generator & rng, cache_index & ci,
                   std::pair<double,double> energy_window, double width_min, double weight_min) :
  INIT_EU_BASE(data,rng,ci,energy_window,width_min,weight_min),
  energy_window(energy_window), width_min(width_min)
 {}

 double attempt() {

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "* Proposing update_change_width (B)" << std::endl;
#endif

  int t = eu::rng(eu::data.temp_conf.size());
  auto const& rect = eu::data.temp_conf[t];

  double dw_min = width_min - rect.width;
  double dw_max = std::min(2*(rect.center - energy_window.first), 2*(energy_window.second - rect.center)) - rect.width;
  double dw = eu::generate_parameter_change(dw_min, dw_max);
  if(dw == 0) return 0;

  eu::update[eu::full].change_rectangle(t, {rect.center, rect.width + dw, rect.height * (1 - dw / (rect.width + dw)), eu::ci});
  eu::update[eu::half].change_rectangle(t, {rect.center, rect.width + dw/2, rect.height * (1 - (dw/2) / (rect.width + dw/2)), eu::ci});

  auto dw_opt = eu::optimize_parameter_change(dw, dw_min, dw_max);

#ifdef EXT_DEBUG
  std::cerr << "Selected rectangle: " << rect << " [" << t << "]" << std::endl;
  std::cerr << "dw_min = " << dw_min << ", dw_max = " << dw_max
            << ", dw = " << dw << ", dw_opt = " << dw_opt.second << std::endl;
#endif

  if(dw_opt.first) {
   dw  = dw_opt.second;
   eu::update[eu::opt].change_rectangle(t, {rect.center, rect.width + dw_opt.second, rect.height * (1 - dw / (rect.width + dw)), eu::ci});
  }

  eu::select_parameter_change(dw_opt.first);

#ifdef EXT_DEBUG
  std::cerr << "selected_parameter_change = " << eu::selected_parameter_change << std::endl;
#endif

  return eu::transition_probability();
 }

};

//-------------------------------------------------------------------------------------

// Change of weight of two rectangles (elementary update C)
template<typename KernelType> class update_change_weight2 : public elementary_update<KernelType> {

 double weight_min;

 using eu = elementary_update<KernelType>;

public:

 update_change_weight2(mc_data<KernelType> & data, random_generator & rng, cache_index & ci,
                       std::pair<double,double> const& energy_window, double width_min, double weight_min) :
  INIT_EU_BASE(data,rng,ci,energy_window,width_min,weight_min),
  weight_min(weight_min)
 {}

 double attempt() {

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "* Proposing update_change_weight2 (C)" << std::endl;
#endif

  int size = eu::data.temp_conf.size();
  if(size < 2) {
#ifdef EXT_DEBUG
   std::cerr << "Not enough rectangles to change" << std::endl;
#endif
   return 0;
  }

  int t1 = eu::rng(size), t2;
  while((t2 = eu::rng(size)) == t1);
  auto const& rect1 = eu::data.temp_conf[t1];
  auto const& rect2 = eu::data.temp_conf[t2];

  double dh1_min = weight_min / rect1.width - rect1.height;
  double dh1_max = (rect2.norm() - weight_min) / rect1.width;
  double dh1 = eu::generate_parameter_change(dh1_min, dh1_max);
  if(dh1 == 0) return 0;

  eu::update[eu::full].change_rectangle(t1, {rect1.center, rect1.width, rect1.height + dh1, eu::ci});
  eu::update[eu::full].change_rectangle(t2, {rect2.center, rect2.width, rect2.height - dh1 * (rect1.width / rect2.width), eu::ci});

  eu::update[eu::half].change_rectangle(t1, {rect1.center, rect1.width, rect1.height + dh1/2, eu::ci});
  eu::update[eu::half].change_rectangle(t2, {rect2.center, rect2.width, rect2.height - (dh1/2) * (rect1.width / rect2.width), eu::ci});

  auto dh1_opt = eu::optimize_parameter_change(dh1, dh1_min, dh1_max);

#ifdef EXT_DEBUG
  std::cerr << "Selected rectangles: " << rect1 << " [" << t1 << "]"
            << " and " << rect2 << " [" << t2 << "]" << std::endl;
  std::cerr << "dh1_min = " << dh1_min << ", dh1_max = " << dh1_max
            << ", dh1 = " << dh1 << ", dh1_opt = " << dh1_opt.second << std::endl;
#endif

  if(dh1_opt.first) {
   dh1  = dh1_opt.second;
   eu::update[eu::opt].change_rectangle(t1, {rect1.center, rect1.width, rect1.height + dh1, eu::ci});
   eu::update[eu::opt].change_rectangle(t2, {rect2.center, rect2.width, rect2.height - dh1 * (rect1.width / rect2.width), eu::ci});
  }

  eu::select_parameter_change(dh1_opt.first);

#ifdef EXT_DEBUG
  std::cerr << "selected_parameter_change = " << eu::selected_parameter_change << std::endl;
#endif

  return eu::transition_probability();
 }

};
}
