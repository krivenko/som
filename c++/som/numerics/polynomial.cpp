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

#include <algorithm>
#include <utility>

#include <triqs/utility/numeric_ops.hpp>

#include "polynomial.hpp"

namespace som {

using nda::range;

////////////////
// polynomial //
////////////////

template <typename CoeffType, bool AutoLowerDegree>
polynomial<CoeffType, AutoLowerDegree>::polynomial(coeffs_type coeffs)
   : coeffs_(std::move(coeffs)) {
  if constexpr(AutoLowerDegree) try_lower_degree();
}

template <typename CoeffType, bool AutoLowerDegree>
polynomial<CoeffType, AutoLowerDegree>::polynomial(
    std::initializer_list<CoeffType> const& coeffs)
   : coeffs_(coeffs) {
  if constexpr(AutoLowerDegree) try_lower_degree();
}

template <typename CoeffType, bool AutoLowerDegree>
polynomial<CoeffType, AutoLowerDegree>&
polynomial<CoeffType, AutoLowerDegree>::operator+=(polynomial const& p) {
  int s1 = coeffs_.size(), s2 = p.coeffs_.size();
  if(s1 > s2)
    coeffs_(range(s2)) += p.coeffs_;
  else if(s1 < s2) {
    coeffs_type tmp(p.coeffs_);
    tmp(range(s1)) += coeffs_;
    std::swap(coeffs_, tmp);
  } else // coeffs_.size() == p.coeffs_.size()
    coeffs_ += p.coeffs_;
  if constexpr(AutoLowerDegree) try_lower_degree();
  return *this;
}

template <typename CoeffType, bool AutoLowerDegree>
polynomial<CoeffType, AutoLowerDegree>&
polynomial<CoeffType, AutoLowerDegree>::operator-=(polynomial const& p) {
  int s1 = coeffs_.size(), s2 = p.coeffs_.size();
  if(s1 > s2)
    coeffs_(range(s2)) -= p.coeffs_;
  else if(s1 < s2) {
    coeffs_type tmp(-p.coeffs_);
    tmp(range(s1)) += coeffs_;
    std::swap(coeffs_, tmp);
  } else // coeffs_.size() == p.coeffs_.size()
    coeffs_ -= p.coeffs_;
  if constexpr(AutoLowerDegree) try_lower_degree();
  return *this;
}

template <typename CoeffType, bool AutoLowerDegree>
polynomial<CoeffType, AutoLowerDegree>&
polynomial<CoeffType, AutoLowerDegree>::operator*=(polynomial const& p) {
  int s1 = coeffs_.size(), s2 = p.coeffs_.size();
  if(s1 == 0 || s2 == 0) {
    coeffs_.resize(0);
    return *this;
  }

  coeffs_type tmp(coeffs_);
  coeffs_.resize(s1 + s2 - 1);
  coeffs_() = CoeffType{};
  for(int n = 0; n < s1; ++n)
    for(int m = 0; m < s2; ++m) { coeffs_(n + m) += tmp(n) * p.coeffs_(m); }
  if constexpr(AutoLowerDegree) try_lower_degree();
  return *this;
}

template <typename CoeffType, bool AutoLowerDegree>
polynomial<CoeffType, AutoLowerDegree>&
polynomial<CoeffType, AutoLowerDegree>::operator*=(CoeffType const& a) {
  int s = coeffs_.size();
  if(s == 0) {
    coeffs_.resize(0);
    return *this;
  }
  coeffs_ *= a;
  if constexpr(AutoLowerDegree) try_lower_degree();
  return *this;
}

template <typename CoeffType, bool AutoLowerDegree>
bool polynomial<CoeffType, AutoLowerDegree>::try_lower_degree(
    double tolerance) {
  using triqs::utility::is_zero;

  int n = coeffs_.size() - 1;
  for(; n >= 0 && is_zero(coeffs_[n], tolerance); --n)
    ;
  if(n != coeffs_.size() - 1) {
    coeffs_ = make_regular(coeffs_(range(n + 1)));
    return true;
  } else
    return false;
}

//////////////////
// derivative() //
//////////////////

template <typename CoeffType, bool AutoLowerDegree>
polynomial<CoeffType, AutoLowerDegree>
derivative(polynomial<CoeffType, AutoLowerDegree> const& p) {
  int s = p.size();
  if(s == 0) return {};

  auto const& p_coeffs = p.coeffs();
  polynomial<CoeffType, AutoLowerDegree> res;
  auto& res_coeffs = res.coeffs();
  res_coeffs.resize(s - 1);
  for(int n = 0; n < s - 1; ++n)
    res_coeffs(n) = p_coeffs(n + 1) * CoeffType(n + 1);
  return res;
}

//////////////////////
// antiderivative() //
//////////////////////

template <typename CoeffType, bool AutoLowerDegree>
polynomial<CoeffType, AutoLowerDegree>
antiderivative(polynomial<CoeffType, AutoLowerDegree> const& p) {
  int s = p.size();
  if(s == 0) return {};

  auto const& p_coeffs = p.coeffs();
  polynomial<CoeffType, AutoLowerDegree> res;
  auto& res_coeffs = res.coeffs();
  res_coeffs.resize(s + 1);
  res_coeffs(0) = CoeffType{};
  for(int n = 1; n < s + 1; ++n) res_coeffs(n) = p_coeffs(n - 1) / CoeffType(n);
  return res;
}

/////////////////////////////
// Explicit instantiations //
/////////////////////////////

template class polynomial<double, true>;
template class polynomial<double, false>;
template class polynomial<std::complex<double>, true>;
template class polynomial<std::complex<double>, false>;

#define INSTANTIATE_F(F, COEFF_TYPE, AUTO_LOWER_DEGREE)                        \
  template polynomial<COEFF_TYPE, AUTO_LOWER_DEGREE> F(                        \
      polynomial<COEFF_TYPE, AUTO_LOWER_DEGREE> const&)

INSTANTIATE_F(derivative, double, true);
INSTANTIATE_F(derivative, double, false);
INSTANTIATE_F(derivative, std::complex<double>, true);
INSTANTIATE_F(derivative, std::complex<double>, false);

INSTANTIATE_F(antiderivative, double, true);
INSTANTIATE_F(antiderivative, double, false);
INSTANTIATE_F(antiderivative, std::complex<double>, true);
INSTANTIATE_F(antiderivative, std::complex<double>, false);

#undef INSTANTIATE_F

} // namespace som
