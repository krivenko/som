/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2026 Igor Krivenko
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

#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>

#include "mesh_traits.hpp"

// FermionGf kernels
#include <som/kernels/fermiongf_imfreq.hpp>
#include <som/kernels/fermiongf_imtime.hpp>
#include <som/kernels/fermiongf_legendre.hpp>

// FermionGfSymm kernels
#include <som/kernels/fermiongfsymm_imfreq.hpp>
#include <som/kernels/fermiongfsymm_imtime.hpp>
#include <som/kernels/fermiongfsymm_legendre.hpp>

// BosonCorr kernels
#include <som/kernels/bosoncorr_imfreq.hpp>
#include <som/kernels/bosoncorr_imtime.hpp>
#include <som/kernels/bosoncorr_legendre.hpp>

// BosonAutoCorr kernels
#include <som/kernels/bosonautocorr_imfreq.hpp>
#include <som/kernels/bosonautocorr_imtime.hpp>
#include <som/kernels/bosonautocorr_legendre.hpp>

// ZeroTemp kernels
#include <som/kernels/zerotemp_imfreq.hpp>
#include <som/kernels/zerotemp_imtime.hpp>
#include <som/kernels/zerotemp_legendre.hpp>

// Declare an external class template 'NAME' by substituting all
// kernel types into its template argument.
#define EXTERN_TEMPLATE_CLASS_FOR_EACH_KERNEL_IMPL(_, ARGS)                    \
  extern template class BOOST_PP_SEQ_ELEM(                                     \
      0, ARGS)<kernel<BOOST_PP_SEQ_ELEM(1, ARGS),                              \
                      triqs::mesh::BOOST_PP_SEQ_ELEM(2, ARGS)>>;
#define EXTERN_TEMPLATE_CLASS_FOR_EACH_KERNEL(NAME)                            \
  BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXTERN_TEMPLATE_CLASS_FOR_EACH_KERNEL_IMPL,    \
                                ((NAME))(ALL_OBSERVABLES)(ALL_INPUT_MESHES))


// Explicitly instantiate a class template 'NAME' by substituting all
// kernel types into its template argument.
#define INSTANTIATE_CLASS_FOR_EACH_KERNEL_IMPL(_, ARGS)                        \
  template class BOOST_PP_SEQ_ELEM(                                            \
      0, ARGS)<kernel<BOOST_PP_SEQ_ELEM(1, ARGS),                              \
                      triqs::mesh::BOOST_PP_SEQ_ELEM(2, ARGS)>>;
#define INSTANTIATE_CLASS_FOR_EACH_KERNEL(NAME)                                \
  BOOST_PP_SEQ_FOR_EACH_PRODUCT(INSTANTIATE_CLASS_FOR_EACH_KERNEL_IMPL,        \
                                ((NAME))(ALL_OBSERVABLES)(ALL_INPUT_MESHES))
