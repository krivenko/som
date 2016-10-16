chi_filename = "chi.h5"
chi_ed_filename = "chi.ed.h5"

gf_struct = {'up':[0], 'dn':[0]}

# Inverse temperature
beta = 40.0

# Impurity parameters
U = 2.0         # Coulomb repulsion
ed = -U*0.6     # Local impurity level
h = 0.1         # Local magnetic field

# Bath parameters
eps = [-1.5,-1.0,-0.5,0.5,1.0,1.5]
V   = [0.3,0.5,0.7,0.8,0.6,0.4]

n_iw = 1024
n_tau = 1001
n_l = 50
