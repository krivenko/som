from h5 import HDFArchive
from triqs.gf import *
from triqs.gf.descriptors import *
from triqs.plot.mpl_interface import oplot
from matplotlib import pyplot as plt
import numpy as np

arch = HDFArchive('consistent_constraints.h5','r')

D = arch['input']['D']
energy_window = arch['accumulate_params']['energy_window']

g_w_som = arch['som_output']['g_w']
chi_som = arch['som_output']['chi2'][0]

g_w_socc = arch['socc_output']['g_w']
chi_socc = arch['socc_output']['chi2'][0]

g_w_ref = g_w_som.copy()
g_w_ref << SemiCircular(D)

w = np.array([float(w) for w in g_w_ref.mesh])
default_model = arch['socc_output']['params']['default_model']

# Plot referenced spectral function
oplot(g_w_ref, mode = 'S', label = "Reference")

# Plot default model
plt.plot(w, default_model, label = "SOCC default model")

# Plot SOM spectral function
oplot(g_w_som, mode = 'S', label = f"SOM, $\\chi^2=%f$" % chi_som)

# Plot SOCC spectral function
oplot(g_w_socc, mode = 'S', label = f"SOCC, $\\chi^2=%f$" % chi_socc)

plt.xlim(energy_window)
plt.xlabel(r"$\omega$")
plt.ylabel(r"$A(\omega)$")

plt.savefig('consistent_constraints.pdf')
