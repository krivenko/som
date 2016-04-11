#pragma once

namespace som {

// All the arguments of the run() function
struct run_parameters_t {

 /// Estimated lower and upper bounds of the spectrum
 /// Negative values of the lower bound will be reset to 0 for susceptibilities and conductivity
 std::pair<double,double> energy_window;

 /// Seed for random number generator
 /// default: 34788 + 928374 * MPI.rank
 int random_seed = 34788 + 928374 * triqs::mpi::communicator().rank();

 /// Name of random number generator
 /// type: str
 std::string random_name = "";

 /// Maximum number of rectangles to represent spectra (:math:`K_{max}`), should be below 70
 unsigned int max_rects = 60;

 /// Minimal width of a rectangle, in units of the energy window width
 double min_rect_width = 1e-3;

 /// Minimal weight of a rectangle, in units of the requested norm of a solution
 double min_rect_weight = 1e-3;

 /// Number of elementary updates per global update (:math:`T`)
 int n_elementary_updates = 1000;

 /// Maximal parameter of the power-law distribution function for the Metropolis algorithm
 double distrib_d_max = 2;

 /// Proposal probability parameter :math:`\gamma`
 double gamma = 2;

 /// Number of global updates (:math:`F`)
 unsigned int n_global_updates = 100;

 /// Adjust the number of global updates automatically
 /// If `true`, use n_global_updates as a starting value
 bool adjust_ngu = true;

 /// Number of particular solutions used to adjust :math:`F`
 unsigned int adjust_ngu_n_solutions = 10;

 /// Limiting value of :math:`\kappa` used to adjust :math:`F`
 double adjust_ngu_kappa = 0.25;

 /// Verbosity level
 /// default: 3 on MPI rank 0, 0 otherwise.
 int verbosity = ((triqs::mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes

 run_parameters_t() {}

 run_parameters_t(std::pair<double,double> energy_window) : energy_window(energy_window) {}

};

}
