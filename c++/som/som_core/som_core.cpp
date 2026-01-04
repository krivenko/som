/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2026 Igor Krivenko
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
#include <type_traits>

#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>

#include "common.hxx"
#include "som_core.hpp"

#include <triqs/utility/exceptions.hpp>

namespace som {

using std::to_string;
using triqs::stat::histogram;

//////////////////////
// som_core::data_t //
//////////////////////

template <typename KernelType>
objective_function<KernelType>
som_core::data_t::make_objf(KernelType const& kernel) const {
  using mesh_t = typename KernelType::mesh_type;
  return std::visit(
      [&kernel, this](auto const& rhs_arg, auto const& errors_arg)
          -> objective_function<KernelType> {
        using rhs_arg_t = std::decay_t<decltype(rhs_arg)>;
        using errors_arg_t = std::decay_t<decltype(errors_arg)>;
        if constexpr(std::is_same_v<rhs_arg_t, input_data_t<mesh_t>> &&
                     std::is_same_v<errors_arg_t, input_data_t<mesh_t>>)
          // Error bars
          return {kernel, rhs_arg, errors_arg};
        else if constexpr(std::is_same_v<rhs_arg_t, input_data_t<mesh_t>> &&
                          std::is_same_v<errors_arg_t, cov_matrix_t<mesh_t>>)
          // Covariance matrix
          return {kernel, rhs_arg, errors_arg, filtering_level};
        else
          TRIQS_RUNTIME_ERROR << "Inconsistent types of input data for "
                                 "objective_function construction";
      },
      rhs,
      errors);
}

#define INSTANTIATE_MAKE_OBJF_IMPL(_, ARGS)                                    \
  template class objective_function<                                           \
      kernel<BOOST_PP_SEQ_ELEM(0, ARGS),                                       \
             triqs::mesh::BOOST_PP_SEQ_ELEM(1, ARGS)>>                         \
  som_core::data_t::make_objf(                                                 \
      kernel<BOOST_PP_SEQ_ELEM(0, ARGS),                                       \
             triqs::mesh::BOOST_PP_SEQ_ELEM(1, ARGS)> const&) const;
// cppcheck-suppress unknownMacro
BOOST_PP_SEQ_FOR_EACH_PRODUCT(INSTANTIATE_MAKE_OBJF_IMPL,
                              (ALL_OBSERVABLES)(ALL_INPUT_MESHES))

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
  std::transform(data.begin(),
                 data.end(),
                 std::back_inserter(conf),
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
  std::transform(data.begin(), data.end(), objf.begin(), [](auto const& d) {
    return d.objf_final;
  });
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
    std::transform(data.begin(),
                   data.end(),
                   std::back_inserter(histograms),
                   [](auto const& d) { return *d.histogram; });
    return {histograms};
  } else
    return {};
}

std::vector<double> som_core::get_objf_min() const {
  std::vector<double> objf_min(data.size());
  std::transform(data.begin(), data.end(), objf_min.begin(), [](auto const& d) {
    return d.objf_min;
  });
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
