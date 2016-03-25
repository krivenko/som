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
#include <stack>

namespace som {

class unique_id;

class id_pool {

 std::stack<int> spare_ids;
 std::vector<int> refcounts;
 const int bucket_size;
 int n_buckets;

 void extend_pool();

public:

 id_pool(int size);
 unique_id get();
 inline void release(int id);
 inline void incref(int id);
 inline void decref(int id);
};

extern id_pool global_id_pool;

class unique_id {

 int id;

 friend class id_pool;

 unique_id(int id);

public:

 ~unique_id();
 unique_id(unique_id const& other);
 unique_id(unique_id && other);
 unique_id & operator=(unique_id const& other);
 unique_id & operator=(unique_id && other);

 operator int();
};

}
