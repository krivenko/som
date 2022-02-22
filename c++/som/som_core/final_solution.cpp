/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#include <cmath>

#include "som_core.hpp"

namespace som {

////////////////////////////////////////
// som_core::compute_final_solution() //
////////////////////////////////////////

void som_core::compute_final_solution(double good_chi_rel,
                                      double good_chi_abs) {

  for(int n = 0; n < data.size(); ++n) {
    auto& d = data[n];

    configuration sol_sum(ci);

    // Rank-local stage of summation
    int n_good_solutions = 0;
    double chi_min = std::sqrt(d.objf_min);
    for(auto const& s : d.particular_solutions) {
      // Pick only good solutions
      double chi = std::sqrt(s.second);
      if(chi / chi_min <= good_chi_rel && chi <= good_chi_abs) {
        sol_sum += s.first;
        ++n_good_solutions;
      }
    }

    // Sum over all processes
    n_good_solutions = mpi::all_reduce(n_good_solutions);
    d.final_solution = mpi::all_reduce(sol_sum);
    d.final_solution *= 1.0 / double(n_good_solutions);
  }
}

/////////////////////
// som_core::run() //
/////////////////////

void som_core::run(accumulate_parameters_t const& p) {
  accumulate(p);
  compute_final_solution(p.adjust_l_good_chi);
}

} // namespace som
