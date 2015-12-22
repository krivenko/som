#include <triqs/test_tools/arrays.hpp>

#include "configuration.hpp"

#define EXPECT_PRINT(X, Y) {std::stringstream ss; ss << Y; EXPECT_EQ(X,ss.str()); }
#define ASSERT_PRINT(X, Y) {std::stringstream ss; ss << Y; ASSERT_EQ(X,ss.str()); }

using namespace som;

configuration conf1 = {{-2.0,2.6,0.3},
                       {1.3,2.6,0.6},
                       {2.0,2.6,0.7}};
std::string conf1_str =
"(c:-2, w:2.6, h:0.3),(c:1.3, w:2.6, h:0.6),(c:2, w:2.6, h:0.7)";
configuration conf2 = {{-3.0,2.6,0.3},
                       {-1.6,2.4,0.1},
                       {2.8,2.8,0.55}};
std::string conf2_str =
"(c:-3, w:2.6, h:0.3),(c:-1.6, w:2.4, h:0.1),(c:2.8, w:2.8, h:0.55)";

TEST(configuration,Print) {
 ASSERT_PRINT(conf1_str,conf1);
 ASSERT_PRINT(conf2_str,conf2);
}

TEST(configuration,Arithmetics) {
 EXPECT_PRINT(conf1_str + "," + conf2_str, conf1 + conf2);

 EXPECT_THROW(-2.0*conf1,triqs::runtime_error);
 EXPECT_THROW(conf1*(-2.0),triqs::runtime_error);

 EXPECT_PRINT("(c:-2, w:2.6, h:0.6),(c:1.3, w:2.6, h:1.2),(c:2, w:2.6, h:1.4)", 2.0*conf1);
 EXPECT_PRINT("(c:-3, w:2.6, h:0.6),(c:-1.6, w:2.4, h:0.2),(c:2.8, w:2.8, h:1.1)", conf2*2.0);

 conf1 += conf2;
 EXPECT_PRINT(conf1_str + "," + conf2_str, conf1);
 conf2 *= 3.0;
 EXPECT_PRINT("(c:-3, w:2.6, h:0.9),(c:-1.6, w:2.4, h:0.3),(c:2.8, w:2.8, h:1.65)", conf2);
}

TEST(configuration,normalize) {
  conf2.normalize(2.56);
  ASSERT_PRINT(conf2_str,conf2);
}

MAKE_MAIN;
