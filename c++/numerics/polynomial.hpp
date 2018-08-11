/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2018 by I. Krivenko
 *
 * SOM is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * SOM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * SOM. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <iostream>
#include <algorithm>
#include <boost/operators.hpp>
#include <triqs/arrays/vector.hpp>
#include <triqs/utility/numeric_ops.hpp>

namespace som {

using triqs::arrays::vector;
using triqs::arrays::range;

// AutoLowerDegree: call try_lower_degree() after every non-const operation
template<typename CoeffType = double, bool AutoLowerDegree = false>
class polynomial : boost::operators<polynomial<CoeffType,AutoLowerDegree>> {

 // Coefficients of the polynomiall; the constant term goes first
 vector<CoeffType> coeffs;

 constexpr static double auto_lower_degree_tolerance = 1e-16;

public:

 polynomial(vector<CoeffType> const& coeffs = {}) : coeffs(coeffs) {
  if(AutoLowerDegree) try_lower_degree();
 }

 polynomial(std::initializer_list<CoeffType> const& coeffs) : coeffs(coeffs) {
  if(AutoLowerDegree) try_lower_degree();
 }

 // Evaluation using Horner's rule
 template<typename ArgType>
 decltype(CoeffType{}*ArgType{}) operator()(ArgType x) const {
  decltype(CoeffType{}*ArgType{}) b = 0;
  for(int i = coeffs.size()-1; i>=0; --i) b = b*x + coeffs[i];
  return b;
 }

 int size() const { return coeffs.size(); }
 int degree() const { return size() - 1; }

 CoeffType  operator[](int n) const { return coeffs[n]; }
 CoeffType& operator[](int n) { return coeffs[n]; }

 bool operator==(const polynomial& p) const {
  if(coeffs.size() == 0) return p.coeffs.size() == 0;
  return coeffs == p.coeffs;
 }

 polynomial & operator+=(polynomial const& p) {
  int s1 = coeffs.size(), s2 = p.coeffs.size();
  if(s1 > s2)
   coeffs(range(s2)) += p.coeffs;
  else if(s1 < s2) {
   vector<CoeffType> tmp(p.coeffs);
   tmp(range(s1)) += coeffs;
   swap(coeffs, tmp);
  } else // coeffs.size() == p.coeffs.size()
   coeffs += p.coeffs;
  if(AutoLowerDegree) try_lower_degree();
  return *this;
 }

 polynomial & operator-=(polynomial const& p) {
  int s1 = coeffs.size(), s2 = p.coeffs.size();
  if(s1 > s2)
   coeffs(range(s2)) -= p.coeffs;
  else if(s1 < s2) {
   vector<CoeffType> tmp(-p.coeffs);
   tmp(range(s1)) += coeffs;
   swap(coeffs, tmp);
  } else // coeffs.size() == p.coeffs.size()
   coeffs -= p.coeffs;
  if(AutoLowerDegree) try_lower_degree();
  return *this;
 }

 polynomial & operator*=(polynomial const& p) {
  int s1 = coeffs.size(), s2 = p.coeffs.size();
  if(s1 == 0 || s2 == 0) {
   coeffs.resize(0);
   return *this;
  }

  vector<CoeffType> tmp(coeffs);
  coeffs.resize(s1 + s2 - 1);
  coeffs() = CoeffType{};
  for(int n = 0; n < s1; ++n)
  for(int m = 0; m < s2; ++m) {
   coeffs(n+m) += tmp(n)*p.coeffs(m);
  }
  if(AutoLowerDegree) try_lower_degree();
  return *this;
 }

 polynomial & operator*=(CoeffType const& a) {
  int s = coeffs.size();
  if(s == 0) {
   coeffs.resize(0);
   return *this;
  }
  coeffs *= a;
  if(AutoLowerDegree) try_lower_degree();
  return *this;
 }
 friend polynomial operator*(polynomial const& p, CoeffType const& a) {
  polynomial res(p); res *= a; return res;
 }
 friend polynomial operator*(CoeffType const& a, polynomial const& p) {
  polynomial res(p); res *= a; return res;
 }

 friend polynomial derivative(polynomial const& p) {
  int s = p.size();
  if(s == 0) return {};

  polynomial res;
  res.coeffs.resize(s-1);
  for(int n = 0; n < s-1; ++n)
   res.coeffs(n) = p.coeffs(n+1) * CoeffType(n+1);
  return res;
 }

 friend polynomial antiderivative(polynomial const& p) {
  int s = p.size();
  if(s == 0) return {};

  polynomial res;
  res.coeffs.resize(s+1);
  res.coeffs(0) = CoeffType{};
  for(int n = 1; n < s+1; ++n)
   res.coeffs(n) = p.coeffs(n-1) / CoeffType(n);
  return res;
 }

 bool try_lower_degree(double tolerance = auto_lower_degree_tolerance) {
  using triqs::utility::is_zero;

  int n = coeffs.size() - 1;
  for(; n >= 0 && is_zero(coeffs[n], tolerance); --n);
  if(n != coeffs.size() - 1) {
   coeffs = coeffs(range(n+1));
   return true;
  } else
   return false;
 }

 friend std::ostream & operator<<(std::ostream & os, polynomial const& p) {
  if(p.size() == 0)
   os << "0";
  else
   for(int n = 0; n < p.size(); ++n) {
    os << p.coeffs[n];
    switch(n) {
     case 0: break;
     case 1: os << "*x"; break;
     default: os << "*x^" << n; break;
    }
    if(n < p.size()-1) os << " + ";
   }
  return os;
 }

};

}
