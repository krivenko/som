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

#include <utility>

#include <triqs/utility/numeric_ops.hpp>

#include "solution_worker.hpp"
#include "updates_class1.hpp"
#include "updates_class2.hpp"

namespace som {

using namespace triqs::mc_tools;

///////////////////
// dist_function //
///////////////////

dist_function::dist_function(random_generator& rng, int T, double d_max)
   : rng(rng), T(T), d_max(d_max) {
  reset();
}

double dist_function::operator()(double x) const {
  return std::pow(x, step < T1 ? 1 + d1 : 1 + d2);
}

void dist_function::reset() {
  step = 0;
  T1 = int(rng(1, T));
  d1 = rng(1);
  d2 = rng(1, d_max);
}

void dist_function::operator++() { ++step; }

/////////////////////
// solution_worker //
/////////////////////

template <typename KernelType>
solution_worker<KernelType>::solution_worker(
    objective_function<KernelType> const& objf, double norm, cache_index& ci,
    worker_parameters_t const& params, std::function<bool()> stop_callback, int f)
   : ci(ci)
   , verbose_mc(params.verbosity >= 3)
   , f(f)
   , t(params.t)
   , mc(params.random_name, params.random_seed, verbose_mc ? 3 : 0)
   , stop_callback(std::move(stop_callback))
   , data{objf,
          {{}, ci},
          {{}, ci},
          0,
          0,
          dist_function(mc.get_rng(), params.t, params.distrib_d_max),
          params.gamma}
   , kern(objf.get_kernel())
   , energy_window(params.energy_window)
   , norm(norm)
   , width_min(params.min_rect_width *
               (params.energy_window.second - params.energy_window.first))
   , weight_min(params.min_rect_weight * norm) {

  // Add all elementary updates
  mc.add_move(update_shift<KernelType>(data, mc.get_rng(), ci, energy_window,
                                       width_min, weight_min),
              "update_shift", 1.0);
  mc.add_move(update_change_width<KernelType>(
                  data, mc.get_rng(), ci, energy_window, width_min, weight_min),
              "update_change_width", 1.0);
  mc.add_move(update_change_weight2<KernelType>(
                  data, mc.get_rng(), ci, energy_window, width_min, weight_min),
              "update_change_weight2", 1.0);
  mc.add_move(update_insert<KernelType>(data, mc.get_rng(), ci, energy_window,
                                        width_min, weight_min,
                                        params.max_rects),
              "update_insert", 1.0);
  mc.add_move(update_remove_shift<KernelType>(
                  data, mc.get_rng(), ci, energy_window, width_min, weight_min),
              "update_remove_shift", 1.0);
  mc.add_move(update_split_shift<KernelType>(data, mc.get_rng(), ci,
                                             energy_window, width_min,
                                             weight_min, params.max_rects),
              "update_split_shift", 1.0);
  mc.add_move(update_glue_shift<KernelType>(
                  data, mc.get_rng(), ci, energy_window, width_min, weight_min),
              "update_glue_shift", 1.0);

  // Reset temporary configuration to global_conf after each global update
  mc.set_after_cycle_duty(std::bind(&solution_worker::reset_temp_conf, this));
}

template <typename KernelType>
configuration
solution_worker<KernelType>::operator()(configuration const& init_config) {

  double init_norm = init_config.norm();
  if(!triqs::utility::is_zero(std::abs(init_norm - norm), 1e-14 * norm))
    TRIQS_RUNTIME_ERROR
        << "solution_worker: initial configuration has a wrong norm "
        << init_norm << ", must be " << norm;

  configuration conf(init_config);

#ifdef EXT_DEBUG
  std::cerr << "solution_worker: using initial configuration (size = "
            << conf.size() << ", norm = " << conf.norm() << "):" << std::endl;
  std::cerr << conf << std::endl;
#endif

  run(conf);
  return conf;
}

template <typename KernelType>
configuration solution_worker<KernelType>::operator()(int init_config_size) {

  if(init_config_size * weight_min > norm)
    TRIQS_RUNTIME_ERROR << "solution_worker: requested configuration size "
                        << init_config_size
                        << " is too big (inconsistent with min_rect_weight = "
                        << weight_min << " and norm = " << norm
                        << ")";

#ifdef EXT_DEBUG
  std::cerr
      << "solution_worker: using a randomly generated initial configuration "
      << "(size = " << init_config_size << ", norm = " << norm << ")"
      << std::endl;
#endif

  auto& rng = mc.get_rng();

  configuration conf(ci);
  for(int i = 0; i < init_config_size; ++i) {
    double center = rng(energy_window.first + width_min / 2,
                        energy_window.second - width_min / 2);
    double width =
        rng(width_min, std::min(2 * (center - energy_window.first),
                                2 * (energy_window.second - center)));
    double height = rng(init_config_size * weight_min / width, norm / width);
    conf.insert({center, width, height, ci});
  }
  conf.normalize(norm);

#ifdef EXT_DEBUG
  std::cerr << "Generated configuration (size = " << conf.size()
            << ", norm = " << conf.norm() << "):" << std::endl;
  std::cerr << conf << std::endl;

#ifndef NDEBUG
  for(auto const& r : conf) assert(r.norm() >= weight_min);
#endif
#endif

  run(conf);
  return conf;
}

template <typename KernelType>
void solution_worker<KernelType>::run(configuration& conf) {

#ifdef EXT_DEBUG
  std::cerr << "solution_worker: starting simulation ..." << std::endl;
#endif

  using std::swap;
  swap(data.global_conf, conf);
  kern.cache_swap(data.global_conf, conf);

  data.temp_conf = data.global_conf;
  if(data.global_conf.cache_ptr.deref().valid)
    kern.cache_copy(data.global_conf, data.temp_conf);
  data.temp_objf_value = data.global_objf_value = data.objf(data.global_conf);

  // Start simulation
  data.Z.reset();
  int res_code = mc.accumulate(f, t, stop_callback);

  swap(data.global_conf, conf);
  kern.cache_swap(data.global_conf, conf);

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
  if(verbose_mc) mc.collect_results(MPI_COMM_SELF);

  // Stopped prematurely
  if(res_code) throw(stopped(res_code));

#ifdef EXT_DEBUG
  std::cerr << "solution_worker: simulation ended." << std::endl;
#endif
}

template <typename KernelType>
void solution_worker<KernelType>::reset_temp_conf() {
  data.temp_conf = data.global_conf;
  if(data.global_conf.cache_ptr.deref().valid)
    kern.cache_copy(data.global_conf, data.temp_conf);
  data.temp_objf_value = data.global_objf_value;
  kern.cache_copy(data.global_conf, data.temp_conf);

#ifdef EXT_DEBUG
  std::cerr << "Temporary configuration reset to (D = " << data.temp_objf_value
            << ")" << std::endl;
  std::cerr << data.temp_conf << std::endl;
#endif

  data.Z.reset();
}

INSTANTIATE_CLASS_FOR_EACH_KERNEL(solution_worker)

} // namespace som
