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

 /// Maximum number of rectangles to represent spectra (K_{max}), should be below 70
 unsigned int max_rects = 60;

 /// Minimal width of a rectangle, in units of the energy window width
 double min_rect_width = 1e-3;

 /// Minimal weight of a rectangle, in units of the requested norm of a solution
 double min_rect_weight = 1e-3;

 /// Number of elementary updates per global update (T)
 int n_elementary_updates = 1000;

 /// Maximal parameter of the power-law distribution function for the Metropolis algorithm
 double distrib_d_max = 2;

 /// Proposal probability parameter :math:`\gamma`
 double gamma = 2;

 run_parameters_t() {}

 run_parameters_t(std::pair<double,double> energy_window) : energy_window(energy_window) {}

};

}
