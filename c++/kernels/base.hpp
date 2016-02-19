/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016 by I. Krivenko
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <vector>

#include "../rectangle.hpp"
#include "../configuration.hpp"
#include "../config_update.hpp"

namespace som {

// Kinds of observables
enum observable_kind {FermionGf, Susceptibility, Conductivity};

// Integral kernels
template<observable_kind, typename Mesh> class kernel;

// Base class for all kernels (CRTP)
template<typename Derived, typename ResultType>
class kernel_base {

public:

 kernel_base() = default;

 inline ResultType operator()(rectangle const& rect) const {
  auto const* kern = static_cast<Derived const*>(this);
  ResultType res(kern->mesh.size());
  kern->apply(rect, res);
  return res;
 }

 inline ResultType operator()(configuration const& c) const {
  auto const* kern = static_cast<Derived const*>(this);
  ResultType res(kern->mesh.size());
  res() = 0;
  for(auto const& r : c) res += operator()(r);
  return res;
 }

 inline ResultType operator()(config_update const& cu) const {
  auto const& conf = cu.get_config();
  ResultType res = operator()(conf);
  for(auto const& ch : cu.changed_rects) {
   int index = ch.first;
   if(index == INT_MAX) { // add rectangle
    res += operator()(ch.second);
   } else if(index < 0) { // remove rectangle
    res -= operator()(conf[-index-1]);
   } else { // change rectangle
    res += operator()(ch.second) - operator()(conf[index]);
   }
  }
  return res;
 }

};

}
