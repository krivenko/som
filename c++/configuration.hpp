#pragma once

#include <vector>

namespace som {

struct rectangle {
 double center;
 double width;
 double height;

 rectangle(double center, double width, double height) :
  center(center), width(width), height(height) {}
 double norm() const { return width*height; }
 double operator()(double x) const {
  return (x >= center-width/2 && x <= center+width/2) ? height : 0;
 }
};

struct configuration {
 std::vector<rectangle> rects;

 configuration() = default;
 configuration(std::initializer_list<rectangle> const& l) : rects(l) {}

 std::size_t size() const { return rects.size(); }
 rectangle const& operator[](std::size_t i) { return rects[i]; }
};

}
