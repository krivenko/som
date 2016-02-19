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

#include <exception>
#include <memory>
#include <vector>
#include <stack>

namespace som {

template<typename T>
class pool {

 struct deleter {
  pool* const p;
  const int id;
  deleter() noexcept {}
  deleter(pool<T> * p, int id) noexcept : p(p), id(id) {}
  deleter(deleter const&) noexcept = delete;
  deleter(deleter &&) noexcept = default;
  void operator()(T *) noexcept { p->spare_ids.push(id); }
 };
 friend struct deleter;

 std::vector<T> data;
 std::stack<int> spare_ids;

public:

 using ptr_type = std::unique_ptr<T,deleter>;

 template<typename... Args>
 pool(int size, Args&&... args) : data(size, T(args...)) {
  for(int id = data.size() - 1; id >= 0; --id) spare_ids.push(id);
 }

 ptr_type get() throw(std::bad_alloc) {
  if(spare_ids.size() == 0) throw(std::bad_alloc());
  int id = spare_ids.top();
  spare_ids.pop();
  return ptr_type(&data[id], deleter(this, id));
 }

};

}
