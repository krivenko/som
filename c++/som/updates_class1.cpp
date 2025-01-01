/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2025 Igor Krivenko
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
#include <iostream>
#endif

#include "solution_worker.hpp"
#include "updates_class1.hpp"

namespace som {

//////////////////
// update_shift //
//////////////////

template <typename KernelType> double update_shift<KernelType>::attempt() {
  eu::attempt_cc_update();

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
  std::cerr << "* Proposing update_shift (A)\n";
#endif

  auto& data = eu::get_data();
  auto& rng = eu::get_rng(); // cppcheck-suppress constVariable
  auto& ci = eu::get_ci();   // cppcheck-suppress constVariable

  int t = rng(data.temp_conf.size());
  auto const& rect = data.temp_conf[t];

  double dc_min = energy_window.first - rect.left();
  double dc_max = energy_window.second - rect.right();
  double dc = eu::generate_parameter_change(dc_min, dc_max);
  if(dc == 0) return 0;

  eu::get_update(eu::full).change_rectangle(
      t, {rect.center + dc, rect.width, rect.height, ci});
  eu::get_update(eu::half).change_rectangle(
      t, {rect.center + dc / 2, rect.width, rect.height, ci});

  auto dc_opt = eu::optimize_parameter_change(dc, dc_min, dc_max);

#ifdef EXT_DEBUG
  std::cerr << "Selected rectangle: " << rect << " [" << t << "]\n";
  std::cerr << "dc_min = " << dc_min << ", dc_max = " << dc_max
            << ", dc = " << dc << ", dc_opt = " << dc_opt.second << '\n';
#endif

  if(dc_opt.first) {
    dc = dc_opt.second;
    eu::get_update(eu::opt).change_rectangle(
        t, {rect.center + dc, rect.width, rect.height, ci});
  }

  eu::select_parameter_change(dc_opt.first);

#ifdef EXT_DEBUG
  std::cerr << "selected_parameter_change = "
            << eu::get_selected_parameter_change() << '\n';
#endif

  return eu::transition_probability();
}

/////////////////////////
// update_change_width //
/////////////////////////

// Change of width without change of weight (elementary update B)
template <typename KernelType>
double update_change_width<KernelType>::attempt() {
  eu::attempt_cc_update();

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
  std::cerr << "* Proposing update_change_width (B)\n";
#endif

  auto& data = eu::get_data();
  auto& rng = eu::get_rng(); // cppcheck-suppress constVariable
  auto& ci = eu::get_ci();   // cppcheck-suppress constVariable

  int t = rng(data.temp_conf.size());
  auto const& rect = data.temp_conf[t];

  double dw_min = width_min - rect.width;
  double dw_max = std::min(2 * (rect.center - energy_window.first),
                           2 * (energy_window.second - rect.center)) -
                  rect.width;
  double dw = eu::generate_parameter_change(dw_min, dw_max);
  if(dw == 0) return 0;

  eu::get_update(eu::full).change_rectangle(
      t,
      {rect.center,
       rect.width + dw,
       rect.height * (1 - dw / (rect.width + dw)),
       ci});
  eu::get_update(eu::half).change_rectangle(
      t,
      {rect.center,
       rect.width + dw / 2,
       rect.height * (1 - (dw / 2) / (rect.width + dw / 2)),
       ci});

  auto dw_opt = eu::optimize_parameter_change(dw, dw_min, dw_max);

#ifdef EXT_DEBUG
  std::cerr << "Selected rectangle: " << rect << " [" << t << "]\n";
  std::cerr << "dw_min = " << dw_min << ", dw_max = " << dw_max
            << ", dw = " << dw << ", dw_opt = " << dw_opt.second << '\n';
#endif

  if(dw_opt.first) {
    dw = dw_opt.second;
    eu::get_update(eu::opt).change_rectangle(
        t,
        {rect.center,
         rect.width + dw_opt.second,
         rect.height * (1 - dw / (rect.width + dw)),
         ci});
  }

  eu::select_parameter_change(dw_opt.first);

#ifdef EXT_DEBUG
  std::cerr << "selected_parameter_change = "
            << eu::get_selected_parameter_change() << '\n';
#endif

  return eu::transition_probability();
}

///////////////////////////
// update_change_weight2 //
///////////////////////////

template <typename KernelType>
double update_change_weight2<KernelType>::attempt() {
  eu::attempt_cc_update();

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
  std::cerr << "* Proposing update_change_weight2 (C)\n";
#endif

  auto& data = eu::get_data();
  auto& rng = eu::get_rng(); // cppcheck-suppress constVariable
  auto& ci = eu::get_ci();   // cppcheck-suppress constVariable

  int size = data.temp_conf.size();
  if(size < 2) {
#ifdef EXT_DEBUG
    std::cerr << "Not enough rectangles to change\n";
#endif
    return 0;
  }

  int t1 = rng(size);
  int t2; // NOLINT(cppcoreguidelines-init-variables)
  while((t2 = rng(size)) == t1)
    ;
  auto const& rect1 = data.temp_conf[t1];
  auto const& rect2 = data.temp_conf[t2];

  double dh1_min = weight_min / rect1.width - rect1.height;
  double dh1_max = (rect2.norm() - weight_min) / rect1.width;
  double dh1 = eu::generate_parameter_change(dh1_min, dh1_max);
  if(dh1 == 0) return 0;

  eu::get_update(eu::full).change_rectangle(
      t1, {rect1.center, rect1.width, rect1.height + dh1, ci});
  eu::get_update(eu::full).change_rectangle(
      t2,
      {rect2.center,
       rect2.width,
       rect2.height - dh1 * (rect1.width / rect2.width),
       ci});

  eu::get_update(eu::half).change_rectangle(
      t1, {rect1.center, rect1.width, rect1.height + dh1 / 2, ci});
  eu::get_update(eu::half).change_rectangle(
      t2,
      {rect2.center,
       rect2.width,
       rect2.height - (dh1 / 2) * (rect1.width / rect2.width),
       ci});

  auto dh1_opt = eu::optimize_parameter_change(dh1, dh1_min, dh1_max);

#ifdef EXT_DEBUG
  std::cerr << "Selected rectangles: " << rect1 << " [" << t1 << "]"
            << " and " << rect2 << " [" << t2 << "]\n";
  std::cerr << "dh1_min = " << dh1_min << ", dh1_max = " << dh1_max
            << ", dh1 = " << dh1 << ", dh1_opt = " << dh1_opt.second << '\n';
#endif

  if(dh1_opt.first) {
    dh1 = dh1_opt.second;
    eu::get_update(eu::opt).change_rectangle(
        t1, {rect1.center, rect1.width, rect1.height + dh1, ci});
    eu::get_update(eu::opt).change_rectangle(
        t2,
        {rect2.center,
         rect2.width,
         rect2.height - dh1 * (rect1.width / rect2.width),
         ci});
  }

  eu::select_parameter_change(dh1_opt.first);

#ifdef EXT_DEBUG
  std::cerr << "selected_parameter_change = "
            << eu::get_selected_parameter_change() << '\n';
#endif

  return eu::transition_probability();
}

INSTANTIATE_CLASS_FOR_EACH_KERNEL(update_shift)
INSTANTIATE_CLASS_FOR_EACH_KERNEL(update_change_width)
INSTANTIATE_CLASS_FOR_EACH_KERNEL(update_change_weight2)

} // namespace som
