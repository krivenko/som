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
#include <numeric>

#include <nda/nda.hpp>
#include <nda/h5.hpp>

#include <triqs/utility/exceptions.hpp>

#include "configuration.hpp"

namespace som {

rectangle& configuration::insert(rectangle const& r) {
  rects.push_back(r);
  return rects.back();
}

void configuration::remove(int index) {
  std::swap(rects[index], rects.back());
  rects.pop_back();
}

void configuration::replace(int index, rectangle const& r) { rects[index] = r; }

void configuration::clear() {
  cache_ptr.invalidate_entry();
  rects.clear();
}

configuration::configuration(std::initializer_list<rectangle> const& l)
   : rects(l) {
  rects.reserve(default_max_rects);
}

configuration::configuration(cache_index& ci, int reserved_rects)
   : cache_ptr(ci) {
  rects.reserve(reserved_rects);
}

configuration::configuration(std::initializer_list<rectangle> const& l,
                             cache_index& ci)
   : rects(l), cache_ptr(ci) {
  rects.reserve(default_max_rects);
}

configuration::configuration(configuration const& c)
   : rects(c.rects),
  cache_ptr(c.cache_ptr ?
            cache_entry_ptr(c.cache_ptr.get_ci()) :
            cache_entry_ptr())
{}

configuration::configuration(configuration&& c) noexcept
   : rects(std::move(c.rects)),
  cache_ptr(c.cache_ptr ?
            cache_entry_ptr(c.cache_ptr.get_ci()) :
            cache_entry_ptr())
{}

// cppcheck-suppress operatorEqVarError
configuration& configuration::operator=(configuration const& c) {
  rects = c.rects;
  cache_ptr.invalidate_entry();
  return *this;
}

configuration& configuration::operator=(configuration&& c) noexcept {
  std::swap(rects, c.rects);
  cache_ptr.invalidate_entry();
  return *this;
}

double configuration::norm() const {
  return std::accumulate(
      rects.begin(), rects.end(), double(0),
      [](double n, rectangle const& r) { return n + r.norm(); });
}

double configuration::operator()(double e) const {
  return std::accumulate(
      rects.begin(), rects.end(), double(0),
      [e](double s, rectangle const& r) { return s + r(e); });
}

configuration& configuration::operator+=(configuration const& c) {
  rects.reserve(rects.size() + c.rects.size());
  std::copy(c.rects.begin(), c.rects.end(), std::back_inserter(rects));
  cache_ptr.invalidate_entry();
  return *this;
}
configuration configuration::operator+(configuration const& c) const {
  configuration s(*this);
  s += c;
  return s;
}

configuration& configuration::operator*=(double alpha) {
  std::transform(rects.begin(), rects.end(), rects.begin(),
                 [&](rectangle const& r) { return r * alpha; });
  cache_ptr.invalidate_entry();
  return *this;
}
configuration operator*(configuration const& c, double alpha) {
  configuration p(c);
  p *= alpha;
  return p;
}
configuration operator*(double alpha, configuration const& c) {
  configuration p(c);
  p *= alpha;
  return p;
}

void configuration::normalize(double norm) {
  if(norm <= 0)
    TRIQS_RUNTIME_ERROR << "Configuration must have a positive norm; norm = "
                        << norm << " is requested";
  double old_norm = std::accumulate(
      rects.begin(), rects.end(), .0,
      [](double n, rectangle const& r) { return n + r.norm(); });
  std::transform(rects.begin(), rects.end(), rects.begin(),
                 [&](rectangle const& r) { return r * (norm / old_norm); });
  cache_ptr.invalidate_entry();
}

std::ostream& operator<<(std::ostream& os, configuration const& c) {
  bool print_comma = false;
  for(auto const& r : c.rects) {
    if(print_comma)
      os << ",";
    else
      print_comma = true;
    os << r;
  }
  return os;
}

void h5_write(h5::group gr, std::string const& name,
              configuration const& c) {
  using nda::array;
  array<double, 2> data(c.size(), 3);
  for(int i = 0; i < c.size(); ++i)
    data(i, nda::range()) =
        array<double, 1>{c[i].center, c[i].width, c[i].height};
  h5_write(gr, name, data);
  auto ds = gr.open_dataset(name);
  write_hdf5_format(ds, c);
}

void h5_read(h5::group gr, std::string const& name, configuration& c) {
  nda::array<double, 2> data;
  h5_read(gr, name, data);
  c.clear();
  for(int i = 0; i < first_dim(data); ++i)
    if(c.cache_ptr)
      c.insert({data(i, 0), data(i, 1), data(i, 2), c.cache_ptr.get_ci()});
    else
      c.insert({data(i, 0), data(i, 1), data(i, 2)});
}

} // namespace som
