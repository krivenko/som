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
#include <numeric>

#include "configuration.hpp"

#include "spectral_stats.hpp"

namespace som {

using namespace nda;

// Overlap of two segments [min_1;max_1] and [min_2;max_2]
inline double overlap(double min1, double max1, double min2, double max2) {
  return std::max(.0, std::min(max1, max2) - std::max(min1, min2));
}

template <typename F>
inline double spectral_integral_impl(configuration const& c, F&& f) {
  return std::accumulate(
      c.begin(), c.end(), .0, [&](double s, struct rectangle const& r) {
        return s + f(r);
      });
}

/////////////////////////
// spectral_integral() //
/////////////////////////

double spectral_integral(double z_m,
                         double delta_m,
                         configuration const& c,
                         resolution_function r_func) {
  switch(r_func) {
    case rectangle:
      return spectral_integral_impl(
          c,
          [delta_m, z_min = z_m - delta_m / 2, z_max = z_m + delta_m / 2](
              struct rectangle const& r) {
            return (r.height / delta_m) *
                   overlap(z_min, z_max, r.left(), r.right());
          });
    case lorentzian:
      return spectral_integral_impl(
          c, [z_m, delta_m](struct rectangle const& r) {
            return (r.height / M_PI) *
                   (std::atan((r.right() - z_m) / (delta_m / 2)) -
                    std::atan((r.left() - z_m) / (delta_m / 2)));
          });
    case gaussian:
      return spectral_integral_impl(
          c, [z_m, delta_m](struct rectangle const& r) {
            return (r.height / 2) *
                   (std::erf((r.right() - z_m) / (M_SQRT2 * delta_m / 2)) -
                    std::erf((r.left() - z_m) / (M_SQRT2 * delta_m / 2)));
          });
  }
}

vector<double> spectral_integral(triqs::mesh::refreq const& mesh,
                                 configuration const& c,
                                 resolution_function r_func) {
  vector<double> integrals(mesh.size());
  double de = mesh.delta();
  for(auto e : mesh)
    integrals(e.linear_index()) = spectral_integral(e, de, c, r_func);
  return integrals;
}

vector<double>
spectral_integral(std::vector<std::pair<double, double>> const& intervals,
                  configuration const& c,
                  resolution_function r_func) {
  vector<double> integrals(intervals.size());
  for(int m = 0; m < intervals.size(); ++m) {
    double z_m = (intervals[m].first + intervals[m].second) / 2;
    double delta_m = intervals[m].second - intervals[m].first;
    integrals(m) = spectral_integral(z_m, delta_m, c, r_func);
  }
  return integrals;
}

////////////////////
// spectral_avg() //
////////////////////

vector<double> spectral_avg(som_core const& cont,
                            int i,
                            triqs::mesh::refreq const& mesh,
                            resolution_function r_func) {
  auto const& solutions = cont.get_particular_solutions(i);

  auto avg = vector<double>::zeros({mesh.size()});
  double de = mesh.delta();
  for(auto e : mesh) {
    for(auto const& s : solutions)
      avg(e.linear_index()) += spectral_integral(e, de, s.first, r_func);
  }

  avg = mpi::all_reduce(avg, cont.get_comm());
  int n_solutions = mpi::all_reduce(solutions.size(), cont.get_comm());

  return avg / n_solutions;
}

vector<double>
spectral_avg(som_core const& cont,
             int i,
             std::vector<std::pair<double, double>> const& intervals,
             resolution_function r_func) {
  auto const& solutions = cont.get_particular_solutions(i);

  auto avg = vector<double>::zeros({static_cast<long>(intervals.size())});
  for(int m = 0; m < intervals.size(); ++m) {
    double z_m = (intervals[m].first + intervals[m].second) / 2;
    double delta_m = intervals[m].second - intervals[m].first;
    for(auto const& s : solutions)
      avg(m) += spectral_integral(z_m, delta_m, s.first, r_func);
  }

  avg = mpi::all_reduce(avg, cont.get_comm());
  int n_solutions = mpi::all_reduce(solutions.size(), cont.get_comm());

  return avg / n_solutions;
}

/////////////////////
// spectral_disp() //
/////////////////////

vector<double> spectral_disp(som_core const& cont,
                             int i,
                             triqs::mesh::refreq const& mesh,
                             vector<double> const& avg,
                             resolution_function r_func) {
  auto const& solutions = cont.get_particular_solutions(i);

  auto disp = vector<double>::zeros({mesh.size()});
  double de = mesh.delta();
  for(auto e : mesh) {
    for(auto const& s : solutions) {
      double diff =
          spectral_integral(e, de, s.first, r_func) - avg(e.linear_index());
      disp(e.linear_index()) += diff * diff;
    }
  }

  disp = mpi::all_reduce(disp, cont.get_comm());
  int n_solutions = mpi::all_reduce(solutions.size(), cont.get_comm());

  return disp / n_solutions;
}

vector<double>
spectral_disp(som_core const& cont,
              int i,
              std::vector<std::pair<double, double>> const& intervals,
              vector<double> const& avg,
              resolution_function r_func) {
  auto const& solutions = cont.get_particular_solutions(i);

  auto disp = vector<double>::zeros({static_cast<long>(intervals.size())});
  for(int m = 0; m < intervals.size(); ++m) {
    double z_m = (intervals[m].first + intervals[m].second) / 2;
    double delta_m = intervals[m].second - intervals[m].first;
    for(auto const& s : solutions) {
      double diff = spectral_integral(z_m, delta_m, s.first, r_func) - avg(m);
      disp(m) += diff * diff;
    }
  }

  disp = mpi::all_reduce(disp, cont.get_comm());
  int n_solutions = mpi::all_reduce(solutions.size(), cont.get_comm());

  return disp / n_solutions;
}

/////////////////////
// spectral_corr() //
/////////////////////

// Regular real frequency mesh
matrix<double> spectral_corr(som_core const& cont,
                             int i,
                             triqs::mesh::refreq const& mesh,
                             vector<double> const& avg,
                             resolution_function r_func) {
  auto const& solutions = cont.get_particular_solutions(i);

  auto corr = matrix<double>::zeros({mesh.size(), mesh.size()});
  double de = mesh.delta();
  for(auto e1 : mesh) {
    for(auto e2 : mesh) {
      for(auto const& s : solutions) {
        double diff1 = (spectral_integral(e1, de, s.first, r_func) -
                        avg(e1.linear_index()));
        double diff2 = (spectral_integral(e2, de, s.first, r_func) -
                        avg(e2.linear_index()));
        corr(e1.linear_index(), e2.linear_index()) += diff1 * diff2;
      }
    }
  }

  corr = mpi::all_reduce(corr, cont.get_comm());
  int n_solutions = mpi::all_reduce(solutions.size(), cont.get_comm());

  return corr / n_solutions;
}

matrix<double>
spectral_corr(som_core const& cont,
              int i,
              std::vector<std::pair<double, double>> const& intervals,
              vector<double> const& avg,
              resolution_function r_func) {
  auto const& solutions = cont.get_particular_solutions(i);

  auto corr = matrix<double>::zeros({static_cast<long>(intervals.size()),
                                     static_cast<long>(intervals.size())});
  for(int m1 = 0; m1 < intervals.size(); ++m1) {
    double z_m1 = (intervals[m1].first + intervals[m1].second) / 2;
    double delta_m1 = intervals[m1].second - intervals[m1].first;
    for(int m2 = 0; m2 < intervals.size(); ++m2) {
      double z_m2 = (intervals[m2].first + intervals[m2].second) / 2;
      double delta_m2 = intervals[m2].second - intervals[m2].first;
      for(auto const& s : solutions) {
        double diff1 =
            spectral_integral(z_m1, delta_m1, s.first, r_func) - avg(m1);
        double diff2 =
            spectral_integral(z_m2, delta_m2, s.first, r_func) - avg(m2);
        corr(m1, m2) += diff1 * diff2;
      }
    }
  }

  corr = mpi::all_reduce(corr, cont.get_comm());
  int n_solutions = mpi::all_reduce(solutions.size(), cont.get_comm());

  return corr / n_solutions;
}

} // namespace som
