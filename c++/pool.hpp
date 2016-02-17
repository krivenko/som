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