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

#include "unique_id.hpp"

#define UNIQUE_ID_POOL_BUCKET_SIZE  0xffff

namespace som {

/////////////
// id_pool //
/////////////

// Singleton
id_pool global_id_pool(UNIQUE_ID_POOL_BUCKET_SIZE);

void id_pool::extend_pool() {
 ++n_buckets;
 int id_min = (n_buckets-1)*bucket_size, id_max = n_buckets*bucket_size-1;
 for(int id = id_max; id >= id_min; --id) spare_ids.push(id);
 refcounts.resize(n_buckets*bucket_size, 0);
}


id_pool::id_pool(int size) : bucket_size(size), n_buckets(0) { extend_pool(); }

unique_id id_pool::get() {
 if(spare_ids.size() == 0) extend_pool();
 int id = spare_ids.top();
 spare_ids.pop();
 refcounts[id] = 1;
 return std::move(unique_id(id));
}

inline void id_pool::release(int id) { spare_ids.push(id); }

inline void id_pool::incref(int id) {
 int & refcount = refcounts[id];
 ++refcount;
}
inline void id_pool::decref(int id) {
 int & refcount = refcounts[id];
 --refcount;
 if(refcount == 0) release(id);
}

///////////////
// unique_id //
///////////////

unique_id::unique_id(int id) : id(id) {}

unique_id::~unique_id() { global_id_pool.decref(id); }

unique_id::unique_id(unique_id const& other) : id(other.id) { global_id_pool.incref(id); }

unique_id::unique_id(unique_id && other) : id(other.id) { global_id_pool.incref(id); }

unique_id & unique_id::operator=(unique_id const& other) {
 global_id_pool.decref(id);
 id = other.id;
 global_id_pool.incref(id);
 return *this;
}

unique_id & unique_id::operator=(unique_id && other) {
 std::swap(id,other.id);
 return *this;
}

inline unique_id::operator int() { return id; }

}
