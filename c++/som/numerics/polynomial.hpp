/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2024 Igor Krivenko
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

#include <complex>
#include <initializer_list>
#include <iostream>

#include <boost/operators.hpp>

#include <nda/nda.hpp>

namespace som {

// AutoLowerDegree: call try_lower_degree() after every non-const operation
template <typename CoeffType = double, bool AutoLowerDegree = false>
class polynomial : boost::operators<polynomial<CoeffType, AutoLowerDegree>> {

  // Coefficients of the polynomial; the constant term goes first
  using coeffs_type = nda::vector<CoeffType>;
  coeffs_type coeffs_ = {};

  constexpr static double auto_lower_degree_tolerance = 1e-16;

public:
  using coeff_type = CoeffType;
  static constexpr bool auto_lower_degree = AutoLowerDegree;

  polynomial() = default;
  explicit polynomial(coeffs_type coeffs);
  explicit polynomial(std::initializer_list<CoeffType> const& coeffs);

  // Evaluation using Horner's rule
  template <typename ArgType>
  inline auto operator()(ArgType x) const -> decltype(CoeffType{} * ArgType{}) {
    decltype(CoeffType{} * ArgType{}) b = 0;
    for(int i = coeffs_.size() - 1; i >= 0; --i) b = b * x + coeffs_[i];
    return b;
  }

  [[nodiscard]] coeffs_type const& coeffs() const { return coeffs_; }
  [[nodiscard]] coeffs_type& coeffs() { return coeffs_; }
  [[nodiscard]] int size() const { return coeffs_.size(); }
  [[nodiscard]] int degree() const { return size() - 1; }

  CoeffType operator[](int n) const { return coeffs_[n]; }
  CoeffType& operator[](int n) { return coeffs_[n]; }

  bool operator==(const polynomial& p) const {
    if(coeffs_.size() == 0) return p.coeffs_.size() == 0;
    return coeffs_ == p.coeffs_;
  }

  polynomial& operator+=(polynomial const& p);
  polynomial& operator-=(polynomial const& p);
  polynomial& operator*=(polynomial const& p);
  polynomial& operator*=(CoeffType const& a);
  friend polynomial operator*(polynomial const& p, CoeffType const& a) {
    polynomial res(p);
    res *= a;
    return res;
  }
  friend polynomial operator*(CoeffType const& a, polynomial const& p) {
    polynomial res(p);
    res *= a;
    return res;
  }

  bool try_lower_degree(double tolerance = auto_lower_degree_tolerance);

  friend std::ostream& operator<<(std::ostream& os, polynomial const& p) {
    if(p.size() == 0)
      os << "0";
    else
      for(int n = 0; n < p.size(); ++n) {
        os << p.coeffs_[n];
        switch(n) {
          case 0: break;
          case 1: os << "*x"; break;
          default: os << "*x^" << n; break;
        }
        if(n < p.size() - 1) os << " + ";
      }
    return os;
  }
};

extern template class polynomial<double, true>;
extern template class polynomial<double, false>;
extern template class polynomial<std::complex<double>, true>;
extern template class polynomial<std::complex<double>, false>;

// derivative()

template <typename CoeffType, bool AutoLowerDegree>
polynomial<CoeffType, AutoLowerDegree>
derivative(polynomial<CoeffType, AutoLowerDegree> const& p);

// antiderivative()

template <typename CoeffType, bool AutoLowerDegree>
polynomial<CoeffType, AutoLowerDegree>
antiderivative(polynomial<CoeffType, AutoLowerDegree> const& p);

} // namespace som
