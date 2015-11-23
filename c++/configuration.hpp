#pragma once

#include <vector>

namespace som {

struct rectangle {
 double center;
 double width;
 double height;

 double norm() const { return width*height; }
 double operator()(double x) const {
  return (x >= center-width/2 && x <= center+width/2) ? height : 0;
 }
};

struct configuration {
 std::vector<rectangle> rects;
 
};

}