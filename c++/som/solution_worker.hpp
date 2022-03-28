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

#include <cmath>
#include <exception>
#include <functional>

#include <triqs/mc_tools.hpp>

#include "kernels/all.hpp"

#include "cache_index.hpp"
#include "configuration.hpp"
#include "solution_functionals/objective_function.hpp"
#include "worker_parameters.hpp"

namespace som {

using namespace triqs::mc_tools;

// Exception: worker has been stopped
struct stopped : public std::exception {
  int code;
  explicit stopped(int code) noexcept : code(code) {}
  stopped(stopped const&) noexcept = default;
  stopped(stopped&&) noexcept = default;
  stopped& operator=(stopped const&) noexcept = default;
  stopped& operator=(stopped&&) noexcept = default;
  ~stopped() override = default;

  [[nodiscard]] const char* what() const noexcept override {
    switch(code) {
      case 0: return "Computation has run until the end";
      case 1: return "Computation has run out of time";
      case 2: return "Computation has been stopped by a signal";
    }
    return "";
  }
};

// Distribution function Z(x), see Eq. (46)
class dist_function {

  triqs::mc_tools::random_generator& rng;

  const bool T1_fixed;
  const int T1_max;
  int T1;
  const double d_max;
  double d1, d2;
  int step;

public:
  dist_function(triqs::mc_tools::random_generator& rng, int T, int T1,
                double d_max);

  double operator()(double x) const;
  [[nodiscard]] bool in_stage_a() const;
  void reset();
  void operator++();
};

// Container for MC data
template <typename KernelType> struct mc_data {
  objective_function<KernelType> const& objf; // Objective function
  configuration temp_conf;                    // Temporary configuration
  configuration global_conf;      // Configuration selected by global updates
  double temp_objf_value = NAN;   // \chi^2(temp_conf)
  double global_objf_value = NAN; // \chi^2(global_conf)
  dist_function Z;                // Distribution function Z(x)
  const double gamma = NAN; // Proposal probability parameter :math:`\gamma`
};

template <typename KernelType> class solution_worker {

  cache_index& ci; // LHS cache index

  bool verbose_mc;       // Make MC output statistics and percentage
  int f;                 // Number of global updates
  int t;                 // Number of elementary updates per global update
  mc_generic<double> mc; // Markov chain
  std::function<bool()> stop_callback; // Stop callback function for mc_generic

  // MC data
  mc_data<KernelType> data;

  // Shortcut to the kernel
  KernelType const& kern;

  std::pair<double, double>
      energy_window; // Estimated lower and upper bounds of the spectrum
  double norm;       // Norm of the sought solution
  double width_min;  // Minimal width of a rectangle
  double weight_min; // Minimal weight of a rectangle

public:
  solution_worker(objective_function<KernelType> const& objf,
                  double norm,
                  cache_index& ci,
                  worker_parameters_t const& params,
                  std::function<bool()> stop_callback,
                  int f);

  // Start from a given configuration
  configuration operator()(configuration const& init_config);

  // Start from a random configuration of size init_config_size
  configuration operator()(int init_config_size);

  [[nodiscard]] triqs::mc_tools::random_generator& get_rng() {
    return mc.get_rng();
  }

  [[nodiscard]] double get_objf_value() const { return data.global_objf_value; }

private:
  void run(configuration& conf);

  void reset_temp_conf();
};

EXTERN_TEMPLATE_CLASS_FOR_EACH_KERNEL(solution_worker)

} // namespace som
