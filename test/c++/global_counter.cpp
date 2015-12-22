#include <triqs/test_tools/arrays.hpp>

#include "global_counter.hpp"

using namespace som;

TEST(global_counter,count) {
 triqs::mpi::communicator c;

 global_counter i(c,1);

#define check  c.barrier(); EXPECT_EQ(ref,long(i)); c.barrier();

 long ref = 1; check;

 if(c.rank() == 0) i = 2;
 ref = 2; check;

 if(c.rank() != 0) i = 3;
 ref = (c.size() == 1) ? 2 : 3; check;

 if(c.rank() == 0) i += 5;
 ref += 5; check;

 if(c.rank() != 0) i += 1;
 ref += c.size() - 1; check;

 i -= 3;
 ref -= 3*c.size(); check;

#undef check
}

MAKE_MAIN;
