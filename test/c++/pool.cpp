#include <triqs/test_tools/arrays.hpp>
#include <vector>

#include "pool.hpp"

using namespace som;

struct my_struct {
 static int created_count;
 static int destroyed_count;

 my_struct() { created_count++; }
 my_struct(my_struct const&) { created_count++; }
 ~my_struct() { destroyed_count++; }
};
int my_struct::created_count = 0;
int my_struct::destroyed_count = 0;

TEST(pool,general) {
 pool<my_struct> p(20);
 for(int n = 0; n < 20; ++n) p.get();
 EXPECT_EQ(21, my_struct::created_count);
 EXPECT_EQ(1, my_struct::destroyed_count);
}

TEST(pool,bad_alloc) {
 pool<my_struct> p(10);

 auto make_ms = [&p](int size){
  std::vector<pool<my_struct>::ptr_type> ptrs;
  for(int n = 0; n < size; ++n) ptrs.push_back(p.get());
 };
 EXPECT_THROW(make_ms(20),std::bad_alloc);
}

MAKE_MAIN;
