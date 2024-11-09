/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2024 Igor Krivenko
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
#include "updates_class2.hpp"

namespace som {

///////////////////
// update_insert //
///////////////////

// Adding a new rectangle (elementary update D)
template <typename KernelType> double update_insert<KernelType>::attempt() {
  eu::attempt_cc_update();

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
  std::cerr << "* Proposing update_insert (D)\n";
#endif

  auto& data = eu::get_data();
  auto& rng = eu::get_rng(); // cppcheck-suppress constVariable
  auto& ci = eu::get_ci();   // cppcheck-suppress constVariable

  int size = data.temp_conf.size();
  if(size == max_rects) {
#ifdef EXT_DEBUG
    std::cerr << "Too many rectangles in the temporary configuration\n";
#endif
    return 0;
  }

  int t = rng(size);
  auto const& rect = data.temp_conf[t];
#ifdef EXT_DEBUG
  std::cerr << "Selected rectangle: " << rect << " [" << t << "]\n";
#endif

  if(rect.norm() <= 2 * weight_min) {
#ifdef EXT_DEBUG
    std::cerr << "Selected rectangle is too small\n";
#endif
    return 0;
  }

  double snew_max = rect.norm() - weight_min;

  double cnew = rng(cnew_min, cnew_max);
  double w_max =
      2 * std::min(energy_window.second - cnew, cnew - energy_window.first);
  double wnew = rng(width_min, w_max);

  double snew = eu::generate_parameter_change(weight_min, snew_max);
  double hnew = snew / wnew;

  eu::get_update(eu::full).change_rectangle(
      t, {rect.center, rect.width, rect.height - snew / rect.width, ci});
  eu::get_update(eu::full).add_rectangle({cnew, wnew, hnew, ci});

  if(snew / 2 < weight_min) { // We cannot optimize w.r.t. snew in this case
    eu::set_parameter_change(eu::full);

#ifdef EXT_DEBUG
    std::cerr << "snew_min = " << weight_min << ", snew_max = " << snew_max
              << ", snew = " << snew << '\n';
#endif
  } else {
    eu::get_update(eu::half).change_rectangle(
        t,
        {rect.center, rect.width, rect.height - (snew / 2) / rect.width, ci});
    eu::get_update(eu::half).add_rectangle({cnew, wnew, hnew / 2, ci});

    auto snew_opt = eu::optimize_parameter_change(snew, width_min, snew_max);

#ifdef EXT_DEBUG
    std::cerr << "snew_min = " << weight_min << ", snew_max = " << snew_max
              << ", snew = " << snew << ", snew_opt = " << snew_opt.second
              << '\n';
#endif

    if(snew_opt.first) {
      snew = snew_opt.second;
      hnew = snew / wnew;
      eu::get_update(eu::opt).change_rectangle(
          t, {rect.center, rect.width, rect.height - snew / rect.width, ci});
      eu::get_update(eu::opt).add_rectangle({cnew, wnew, hnew, ci});
    }

    eu::select_parameter_change(snew_opt.first);
  }

#ifdef EXT_DEBUG
  std::cerr << "selected_parameter_change = "
            << eu::get_selected_parameter_change() << '\n';
#endif

  return eu::transition_probability();
}

/////////////////////////
// update_remove_shift //
/////////////////////////

template <typename KernelType>
double update_remove_shift<KernelType>::attempt() {
  eu::attempt_cc_update();

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
  std::cerr << "* Proposing update_remove_shift (E)\n";
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

  double dc2_min = energy_window.first - rect2.left();
  double dc2_max = energy_window.second - rect2.right();
  double dc2 = eu::generate_parameter_change(dc2_min, dc2_max);
  if(dc2 == 0) return 0;

  eu::get_update(eu::full).change_rectangle(
      t2,
      {rect2.center + dc2,
       rect2.width,
       rect2.height + rect1.norm() / rect2.width,
       ci});
  eu::get_update(eu::full).remove_rectangle(t1);
  eu::get_update(eu::half).change_rectangle(
      t2,
      {rect2.center + dc2 / 2,
       rect2.width,
       rect2.height + rect1.norm() / rect2.width,
       ci});
  eu::get_update(eu::half).remove_rectangle(t1);

  auto dc2_opt = eu::optimize_parameter_change(dc2, dc2_min, dc2_max);

#ifdef EXT_DEBUG
  std::cerr << "Selected rectangle for removal: " << rect1 << " [" << t1
            << "]\n";
  std::cerr << "Selected rectangle for shift: " << rect2 << " [" << t2 << "]\n";
  std::cerr << "dc2_min = " << dc2_min << ", dc2_max = " << dc2_max
            << ", dc2 = " << dc2 << ", dc2_opt = " << dc2_opt.second << '\n';
#endif

  if(dc2_opt.first) {
    dc2 = dc2_opt.second;
    eu::get_update(eu::opt).change_rectangle(
        t2,
        {rect2.center + dc2,
         rect2.width,
         rect2.height + rect1.norm() / rect2.width,
         ci});
    eu::get_update(eu::opt).remove_rectangle(t1);
  }

  eu::select_parameter_change(dc2_opt.first);

#ifdef EXT_DEBUG
  std::cerr << "selected_parameter_change = "
            << eu::get_selected_parameter_change() << '\n';
#endif

  return eu::transition_probability();
}

////////////////////////
// update_split_shift //
////////////////////////

template <typename KernelType>
double update_split_shift<KernelType>::attempt() {
  eu::attempt_cc_update();

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
  std::cerr << "* Proposing update_split_shift (F)\n";
#endif

  auto& data = eu::get_data();
  auto& rng = eu::get_rng(); // cppcheck-suppress constVariable
  auto& ci = eu::get_ci();   // cppcheck-suppress constVariable

  int size = data.temp_conf.size();
  if(size == max_rects) {
#ifdef EXT_DEBUG
    std::cerr << "Too many rectangles in the temporary configuration\n";
#endif
    return 0;
  }

  int t = rng(size);
  auto const& rect = data.temp_conf[t];

#ifdef EXT_DEBUG
  std::cerr << "Selected rectangle: " << rect << " [" << t << "]\n";
#endif
  if(rect.width <= 2 * width_min) {
#ifdef EXT_DEBUG
    std::cerr << "Selected rectangle is too narrow\n";
#endif
    return 0;
  }
  if(rect.norm() <= 2 * weight_min) {
#ifdef EXT_DEBUG
    std::cerr << "Selected rectangle is too small\n";
#endif
    return 0;
  }

  double w_min = std::max(width_min, weight_min / rect.height);
  double w_max =
      std::min(rect.width - width_min, rect.width - weight_min / rect.height);
  double w1 = rng(w_min, w_max);
  double w2 = rect.width - w1;
  double c1 = rect.left() + w1 / 2;
  double c2 = rect.right() - w2 / 2;

  if(w1 <= w2) {
    double dc1_min = std::max(energy_window.first + w1 / 2 - c1,
                              (c2 - (energy_window.second - w2 / 2)) * w2 / w1);
    double dc1_max = std::min(energy_window.second - w1 / 2 - c1,
                              (c2 - (energy_window.first + w2 / 2)) * w2 / w1);
    double dc1 = eu::generate_parameter_change(dc1_min, dc1_max);
    if(dc1 == 0) return 0;
    double dc2 = -dc1 * w1 / w2;

    eu::get_update(eu::full).change_rectangle(t,
                                              {c1 + dc1, w1, rect.height, ci});
    eu::get_update(eu::full).add_rectangle({c2 + dc2, w2, rect.height, ci});
    eu::get_update(eu::half).change_rectangle(
        t, {c1 + dc1 / 2, w1, rect.height, ci});
    eu::get_update(eu::half).add_rectangle({c2 + dc2 / 2, w2, rect.height, ci});

    auto dc1_opt = eu::optimize_parameter_change(dc1, dc1_min, dc1_max);

#ifdef EXT_DEBUG
    std::cerr << "dc1_min = " << dc1_min << ", dc1_max = " << dc1_max
              << ", dc1 = " << dc1 << ", dc1_opt = " << dc1_opt.second << '\n';
#endif

    if(dc1_opt.first) {
      dc1 = dc1_opt.second;
      dc2 = -dc1 * w1 / w2;
      eu::get_update(eu::opt).change_rectangle(t,
                                               {c1 + dc1, w1, rect.height, ci});
      eu::get_update(eu::opt).add_rectangle({c2 + dc2, w2, rect.height, ci});
    }

    eu::select_parameter_change(dc1_opt.first);

  } else { // w1 > w2
    double dc2_min = std::max(energy_window.first + w2 / 2 - c2,
                              (c1 - (energy_window.second - w1 / 2)) * w1 / w2);
    double dc2_max = std::min(energy_window.second - w2 / 2 - c2,
                              (c1 - (energy_window.first + w1 / 2)) * w1 / w2);
    double dc2 = eu::generate_parameter_change(dc2_min, dc2_max);
    if(dc2 == 0) return 0;
    double dc1 = -dc2 * w2 / w1;

    eu::get_update(eu::full).change_rectangle(t,
                                              {c1 + dc1, w1, rect.height, ci});
    eu::get_update(eu::full).add_rectangle({c2 + dc2, w2, rect.height, ci});
    eu::get_update(eu::half).change_rectangle(
        t, {c1 + dc1 / 2, w1, rect.height, ci});
    eu::get_update(eu::half).add_rectangle({c2 + dc2 / 2, w2, rect.height, ci});

    auto dc2_opt = eu::optimize_parameter_change(dc2, dc2_min, dc2_max);

#ifdef EXT_DEBUG
    std::cerr << "dc2_min = " << dc2_min << ", dc2_max = " << dc2_max
              << ", dc2 = " << dc2 << ", dc2_opt = " << dc2_opt.second << '\n';
#endif

    if(dc2_opt.first) {
      dc2 = dc2_opt.second;
      dc1 = -dc2 * w2 / w1;
      eu::get_update(eu::opt).change_rectangle(t,
                                               {c1 + dc1, w1, rect.height, ci});
      eu::get_update(eu::opt).add_rectangle({c2 + dc2, w2, rect.height, ci});
    }

    eu::select_parameter_change(dc2_opt.first);
  }

#ifdef EXT_DEBUG
  std::cerr << "selected_parameter_change = "
            << eu::get_selected_parameter_change() << '\n';
#endif

  return eu::transition_probability();
}

///////////////////////
// update_glue_shift //
///////////////////////

template <typename KernelType> double update_glue_shift<KernelType>::attempt() {
  eu::attempt_cc_update();

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
  std::cerr << "* Proposing move_glue_shift (G)\n";
#endif

  auto& data = eu::get_data();
  auto& rng = eu::get_rng(); // cppcheck-suppress constVariable
  auto& ci = eu::get_ci();   // cppcheck-suppress constVariable

  int size = data.temp_conf.size();
  if(size < 2) {
#ifdef EXT_DEBUG
    std::cerr << "Not enough rectangles to glue\n";
#endif
    return 0;
  }

  int t1 = rng(size);
  int t2; // NOLINT(cppcoreguidelines-init-variables)
  while((t2 = rng(size)) == t1)
    ;
  auto const& rect1 = data.temp_conf[t1];
  auto const& rect2 = data.temp_conf[t2];

  double w = (rect1.width + rect2.width) / 2;
  double s1 = rect1.norm(), s2 = rect2.norm();
  double s = s1 + s2;
  double c = (rect1.center * s1 + rect2.center * s2) / s;
  double h = s / w;

  double dc_min = energy_window.first + w / 2 - c;
  double dc_max = energy_window.second - w / 2 - c;
  double dc = eu::generate_parameter_change(dc_min, dc_max);
  if(dc == 0) return 0;

  eu::get_update(eu::full).change_rectangle(t1, {c + dc, w, h, ci});
  eu::get_update(eu::full).remove_rectangle(t2);

#ifdef EXT_DEBUG
  std::cerr << "Selected rectangles: " << rect1 << " [" << t1 << "]"
            << " and " << rect2 << " [" << t2 << "]\n";
#endif

  if(dc_min * dc_max > 0) { // In this case dc/2 can be outside [dc_min; dc_max]
    eu::set_parameter_change(eu::full);

#ifdef EXT_DEBUG
    std::cerr << "dc_min = " << dc_min << ", dc_max = " << dc_max
              << ", dc = " << dc << '\n';
#endif


  } else {
    eu::get_update(eu::half).change_rectangle(t1, {c + dc / 2, w, h, ci});
    eu::get_update(eu::half).remove_rectangle(t2);

    auto dc_opt = eu::optimize_parameter_change(dc, dc_min, dc_max);

#ifdef EXT_DEBUG
    std::cerr << "dc_min = " << dc_min << ", dc_max = " << dc_max
              << ", dc = " << dc << ", dc_opt = " << dc_opt.second << '\n';
#endif

    if(dc_opt.first) {
      dc = dc_opt.second;
      eu::get_update(eu::opt).change_rectangle(t1, {c + dc, w, h, ci});
      eu::get_update(eu::opt).remove_rectangle(t2);
    }

    eu::select_parameter_change(dc_opt.first);
  }

#ifdef EXT_DEBUG
  std::cerr << "selected_parameter_change = "
            << eu::get_selected_parameter_change() << '\n';
#endif

  return eu::transition_probability();
}

INSTANTIATE_CLASS_FOR_EACH_KERNEL(update_insert)
INSTANTIATE_CLASS_FOR_EACH_KERNEL(update_remove_shift)
INSTANTIATE_CLASS_FOR_EACH_KERNEL(update_split_shift)
INSTANTIATE_CLASS_FOR_EACH_KERNEL(update_glue_shift)

} // namespace som
