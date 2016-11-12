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
#pragma once

#include "elementary_update.hpp"

namespace som {

// Adding a new rectangle (elementary update D)
template<typename KernelType> class update_insert : public elementary_update<KernelType> {

 std::pair<double,double> energy_window;
 double width_min;
 double weight_min;
 double cnew_min, cnew_max;
 int max_rects;

 using eu = elementary_update<KernelType>;

public:

 update_insert(mc_data<KernelType> & data, random_generator & rng, cache_index & ci,
                   std::pair<double,double> energy_window, double width_min, double weight_min, int max_rects)  :
  INIT_EU_BASE(data,rng,ci,width_min,weight_min),
  energy_window(energy_window), width_min(width_min), weight_min(weight_min),
  cnew_min(energy_window.first + width_min/2), cnew_max(energy_window.second - width_min/2),
  max_rects(max_rects)
 {}

 double attempt() {

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "* Proposing update_insert (D)" << std::endl;
#endif

  int size = eu::data.temp_conf.size();
  if(size == max_rects) {
#ifdef EXT_DEBUG
   std::cerr << "Too many rectangles in the temporary configuration" << std::endl;
#endif
   return 0;
  }

  int t = eu::rng(size);
  auto const& rect = eu::data.temp_conf[t];
#ifdef EXT_DEBUG
  std::cerr << "Selected rectangle: " << rect << " [" << t << "]" << std::endl;
#endif

  if(rect.norm() <= 2*weight_min) {
#ifdef EXT_DEBUG
   std::cerr << "Selected rectangle is too small" << std::endl;
#endif
   return 0;
  }

  double snew_max = rect.norm() - weight_min;

  double cnew = eu::rng(cnew_min, cnew_max);
  double w_max = 2*std::min(energy_window.second - cnew, cnew - energy_window.first);
  double wnew = eu::rng(width_min, w_max);

  double snew = eu::generate_parameter_change(weight_min, snew_max);
  double hnew = snew / wnew;

  eu::update[eu::full].change_rectangle(t, {rect.center, rect.width, rect.height - snew / rect.width, eu::ci});
  eu::update[eu::full].add_rectangle({cnew, wnew, hnew, eu::ci});

  if(snew / 2 < weight_min) { // We cannot optimize w.r.t. snew in this case
   eu::new_objf_value[eu::full] = eu::data.objf(eu::update[eu::full]);
   eu::selected_parameter_change = eu::full;

#ifdef EXT_DEBUG
  std::cerr << "snew_min = " << weight_min << ", snew_max = " << snew_max
            << ", snew = " << snew << std::endl;
#endif
  } else {
   eu::update[eu::half].change_rectangle(t, {rect.center, rect.width, rect.height - (snew/2) / rect.width, eu::ci});
   eu::update[eu::half].add_rectangle({cnew, wnew, hnew/2, eu::ci});

   auto snew_opt = eu::optimize_parameter_change(snew, width_min, snew_max);

#ifdef EXT_DEBUG
  std::cerr << "snew_min = " << weight_min << ", snew_max = " << snew_max
            << ", snew = " << snew << ", snew_opt = " << snew_opt.second << std::endl;
#endif

   if(snew_opt.first) {
    snew = snew_opt.second;
    hnew = snew / wnew;
    eu::update[eu::opt].change_rectangle(t, {rect.center, rect.width, rect.height - snew / rect.width, eu::ci});
    eu::update[eu::opt].add_rectangle({cnew, wnew, hnew, eu::ci});
   }

   eu::select_parameter_change(snew_opt.first);
  }

#ifdef EXT_DEBUG
  std::cerr << "selected_parameter_change = " << eu::selected_parameter_change << std::endl;
#endif

  return eu::transition_probability();
 }

};


//-------------------------------------------------------------------------------------

// Removing one rectangle and shifting another (elementary update E)
template<typename KernelType> class update_remove_shift : public elementary_update<KernelType> {

 std::pair<double,double> energy_window;

 using eu = elementary_update<KernelType>;

public:

 update_remove_shift(mc_data<KernelType> & data, random_generator & rng, cache_index & ci,
                     std::pair<double,double> energy_window, double width_min, double weight_min) :
  INIT_EU_BASE(data,rng,ci,width_min,weight_min),
  energy_window(energy_window)
 {}

 double attempt() {

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "* Proposing update_remove_shift (E)" << std::endl;
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

  double dc2_min = energy_window.first + rect2.width/2 - rect2.center;
  double dc2_max = energy_window.second - rect2.width/2 - rect2.center;
  double dc2 = eu::generate_parameter_change(dc2_min, dc2_max);
  if(dc2 == 0) return 0;

  eu::update[eu::full].change_rectangle(t2, {rect2.center + dc2, rect2.width, rect2.height + rect1.norm() / rect2.width, eu::ci});
  eu::update[eu::full].remove_rectangle(t1);
  eu::update[eu::half].change_rectangle(t2, {rect2.center + dc2/2, rect2.width, rect2.height + rect1.norm() / rect2.width, eu::ci});
  eu::update[eu::half].remove_rectangle(t1);

  auto dc2_opt = eu::optimize_parameter_change(dc2, dc2_min, dc2_max);

#ifdef EXT_DEBUG
  std::cerr << "Selected rectangle for removal: " << rect1 << " [" << t1 << "]" << std::endl;
  std::cerr << "Selected rectangle for shift: " << rect2 << " [" << t2 << "]" << std::endl;
  std::cerr << "dc2_min = " << dc2_min << ", dc2_max = " << dc2_max
            << ", dc2 = " << dc2 << ", dc2_opt = " << dc2_opt.second << std::endl;
#endif

  if(dc2_opt.first) {
   dc2 = dc2_opt.second;
   eu::update[eu::opt].change_rectangle(t2, {rect2.center + dc2, rect2.width, rect2.height + rect1.norm() / rect2.width, eu::ci});
   eu::update[eu::opt].remove_rectangle(t1);
  }

  eu::select_parameter_change(dc2_opt.first);

#ifdef EXT_DEBUG
  std::cerr << "selected_parameter_change = " << eu::selected_parameter_change << std::endl;
#endif

  return eu::transition_probability();
 }

};

//-------------------------------------------------------------------------------------

// Splitting a rectangle and shifting the pieces (elementary update F)
template<typename KernelType> class update_split_shift : public elementary_update<KernelType> {

 std::pair<double,double> energy_window;
 double width_min;
 double weight_min;
 int max_rects;

 using eu = elementary_update<KernelType>;

public:

 update_split_shift(mc_data<KernelType> & data, random_generator & rng, cache_index & ci,
                   std::pair<double,double> energy_window, double width_min, double weight_min, int max_rects) :
  INIT_EU_BASE(data,rng,ci,width_min,weight_min),
  energy_window(energy_window), width_min(width_min), weight_min(weight_min), max_rects(max_rects)
 {}

 double attempt() {

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "* Proposing update_split_shift (F)" << std::endl;
#endif

  int size = eu::data.temp_conf.size();
  if(size == max_rects) {
#ifdef EXT_DEBUG
   std::cerr << "Too many rectangles in the temporary configuration" << std::endl;
#endif
   return 0;
  }

  int t = eu::rng(size);
  auto const& rect = eu::data.temp_conf[t];

#ifdef EXT_DEBUG
  std::cerr << "Selected rectangle: " << rect << " [" << t << "]" << std::endl;
#endif
  if(rect.width <= 2*width_min) {
#ifdef EXT_DEBUG
   std::cerr << "Selected rectangle is too narrow" << std::endl;
#endif
   return 0;
  }
  if(rect.norm() <= 2*weight_min) {
#ifdef EXT_DEBUG
   std::cerr << "Selected rectangle is too small" << std::endl;
#endif
   return 0;
  }

  double w_min = std::max(width_min, weight_min/rect.height);
  double w_max = std::min(rect.width - width_min, rect.width - weight_min/rect.height);
  double w1 = eu::rng(w_min, w_max);
  double w2 = rect.width - w1;
  double c1 = rect.center - rect.width/2 + w1/2;
  double c2 = rect.center + rect.width/2 - w2/2;

  if(w1 <= w2) {
   double dc1_min = std::max(energy_window.first + w1/2 - c1, (c2 - (energy_window.second - w2/2))*w2/w1);
   double dc1_max = std::min(energy_window.second - w1/2 - c1, (c2 - (energy_window.first + w2/2))*w2/w1);
   double dc1 = eu::generate_parameter_change(dc1_min, dc1_max);
   if(dc1 == 0) return 0;
   double dc2 = -dc1 * w1 / w2;

   eu::update[eu::full].change_rectangle(t, {c1 + dc1, w1, rect.height, eu::ci});
   eu::update[eu::full].add_rectangle({c2 + dc2, w2, rect.height, eu::ci});
   eu::update[eu::half].change_rectangle(t, {c1 + dc1/2, w1, rect.height, eu::ci});
   eu::update[eu::half].add_rectangle({c2 + dc2/2, w2, rect.height, eu::ci});

   auto dc1_opt = eu::optimize_parameter_change(dc1, dc1_min, dc1_max);

#ifdef EXT_DEBUG
   std::cerr << "dc1_min = " << dc1_min << ", dc1_max = " << dc1_max
             << ", dc1 = " << dc1 << ", dc1_opt = " << dc1_opt.second << std::endl;
#endif

   if(dc1_opt.first) {
    dc1 = dc1_opt.second;
    dc2 = -dc1 * w1 / w2;
    eu::update[eu::opt].change_rectangle(t, {c1 + dc1, w1, rect.height, eu::ci});
    eu::update[eu::opt].add_rectangle({c2 + dc2, w2, rect.height, eu::ci});
   }

   eu::select_parameter_change(dc1_opt.first);

  } else { // w1 > w2
   double dc2_min = std::max(energy_window.first + w2/2 - c2, (c1 - (energy_window.second - w1/2))*w1/w2);
   double dc2_max = std::min(energy_window.second - w2/2 - c2, (c1 - (energy_window.first + w1/2))*w1/w2);
   double dc2 = eu::generate_parameter_change(dc2_min, dc2_max);
   if(dc2 == 0) return 0;
   double dc1 = -dc2 * w2 / w1;

   eu::update[eu::full].change_rectangle(t, {c1 + dc1, w1, rect.height, eu::ci});
   eu::update[eu::full].add_rectangle({c2 + dc2, w2, rect.height, eu::ci});
   eu::update[eu::half].change_rectangle(t, {c1 + dc1/2, w1, rect.height, eu::ci});
   eu::update[eu::half].add_rectangle({c2 + dc2/2, w2, rect.height, eu::ci});

   auto dc2_opt = eu::optimize_parameter_change(dc2, dc2_min, dc2_max);

#ifdef EXT_DEBUG
   std::cerr << "dc2_min = " << dc2_min << ", dc2_max = " << dc2_max
             << ", dc2 = " << dc2 << ", dc2_opt = " << dc2_opt.second << std::endl;
#endif

  if(dc2_opt.first) {
    dc2 = dc2_opt.second;
    dc1 = -dc2 * w2 / w1;
    eu::update[eu::opt].change_rectangle(t, {c1 + dc1, w1, rect.height, eu::ci});
    eu::update[eu::opt].add_rectangle({c2 + dc2, w2, rect.height, eu::ci});
   }

   eu::select_parameter_change(dc2_opt.first);
  }

#ifdef EXT_DEBUG
  std::cerr << "selected_parameter_change = " << eu::selected_parameter_change << std::endl;
#endif

  return eu::transition_probability();
 }

};

//-------------------------------------------------------------------------------------

// Glueing rectangles (elementary update G)
template<typename KernelType> class update_glue_shift : public elementary_update<KernelType> {

 std::pair<double,double> energy_window;

 using eu = elementary_update<KernelType>;

public:

 update_glue_shift(mc_data<KernelType> & data, random_generator & rng, cache_index & ci,
                   std::pair<double,double> energy_window, double width_min, double weight_min) :
  INIT_EU_BASE(data,rng,ci,width_min,weight_min),
  energy_window(energy_window)
 {}

 double attempt() {

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "* Proposing move_glue_shift (G)" << std::endl;
#endif

  int size = eu::data.temp_conf.size();
  if(size < 2) {
#ifdef EXT_DEBUG
   std::cerr << "Not enough rectangles to glue" << std::endl;
#endif
   return 0;
  }

  int t1 = eu::rng(size), t2;
  while((t2 = eu::rng(size)) == t1);
  auto const& rect1 = eu::data.temp_conf[t1];
  auto const& rect2 = eu::data.temp_conf[t2];

  double w = (rect1.width + rect2.width) / 2;
  double s1 = rect1.norm(), s2 = rect2.norm();
  double s = s1 + s2;
  double c = (rect1.center * s1 + rect2.center * s2) / s;
  double h = s / w;

  double dc_min = energy_window.first + w/2 - c;
  double dc_max = energy_window.second - w/2 - c;
  double dc = eu::generate_parameter_change(dc_min, dc_max);
  if(dc == 0) return 0;

  eu::update[eu::full].change_rectangle(t1, {c + dc, w, h, eu::ci});
  eu::update[eu::full].remove_rectangle(t2);
  eu::update[eu::half].change_rectangle(t1, {c + dc/2, w, h, eu::ci});
  eu::update[eu::half].remove_rectangle(t2);

  auto dc_opt = eu::optimize_parameter_change(dc, dc_min, dc_max);

#ifdef EXT_DEBUG
  std::cerr << "Selected rectangles: " << rect1 << " [" << t1 << "]"
            << " and " << rect2 << " [" << t2 << "]" << std::endl;
  std::cerr << "dc_min = " << dc_min << ", dc_max = " << dc_max
            << ", dc = " << dc << ", dc_opt = " << dc_opt.second << std::endl;
#endif

  if(dc_opt.first) {
   dc = dc_opt.second;
   eu::update[eu::opt].change_rectangle(t1, {c + dc, w, h, eu::ci});
   eu::update[eu::opt].remove_rectangle(t2);
  }

  eu::select_parameter_change(dc_opt.first);

#ifdef EXT_DEBUG
  std::cerr << "selected_parameter_change = " << eu::selected_parameter_change << std::endl;
#endif

  return eu::transition_probability();
 }

};

}
