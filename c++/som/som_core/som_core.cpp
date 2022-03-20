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

#include <algorithm>
#include <iterator>
#include <limits>

#include "common.hxx"
#include "som_core.hpp"

namespace som {

using std::to_string;
using triqs::stat::histogram;

/////////////////////////
// som_core: Accessors //
/////////////////////////

std::vector<std::pair<configuration, double>> const&
som_core::get_particular_solutions(long i) const {
  if(i >= data.size())
    fatal_error("Matrix element index " + to_string(i) + " out of bounds");
  return data[i].particular_solutions;
}

configuration const& som_core::get_solution(long i) const {
  if(i >= data.size())
    fatal_error("Matrix element index " + to_string(i) + " out of bounds");
  return data[i].final_solution;
}

std::vector<configuration> som_core::get_solutions() const {
  std::vector<configuration> conf;
  conf.reserve(data.size());
  std::transform(data.begin(), data.end(), std::back_inserter(conf),
                 [](auto const& d) { return d.final_solution; });

  return conf;
}

double som_core::get_objf(long i) const {
  if(i >= data.size())
    fatal_error("Matrix element index " + to_string(i) + " out of bounds");
  return data[i].objf_final;
}

std::vector<double> som_core::get_objf() const {
  std::vector<double> objf(data.size());
  std::transform(data.begin(), data.end(), objf.begin(),
                 [](auto const& d) { return d.objf_final; });
  return objf;
}

std::optional<histogram> const& som_core::get_histogram(long i) const {
  if(i >= data.size())
    fatal_error("Matrix element index " + to_string(i) + " out of bounds");
  return data[i].histogram;
}

std::optional<std::vector<histogram>> som_core::get_histograms() const {
  if(data.back().histogram) {
    std::vector<histogram> histograms;
    histograms.reserve(data.size());
    std::transform(data.begin(), data.end(), std::back_inserter(histograms),
    [](auto const& d) { return *d.histogram; });
    return {histograms};
  } else
    return {};
}

std::vector<double> som_core::get_objf_min() const {
  std::vector<double> objf_min(data.size());
  std::transform(data.begin(), data.end(), objf_min.begin(),
                 [](auto const& d) { return d.objf_min; });
  return objf_min;
}

///////////////////////
// som_core::clear() //
///////////////////////

void som_core::clear() {
  for(auto& d : data) {
    d.particular_solutions.clear();
    d.objf_min = HUGE_VAL;
    d.final_solution.clear();
    d.histogram.reset();
  }
}

} // namespace som
